# ── Phylogenetic scaffold pipeline ───────────────────────────────────────────
#
# Overview of the two-stage approach:
#
#   Stage 1 — build_or_read_phylogeny()
#     Uses prepR4pcm to get a tree and produce a name-mapping table.
#     The four-stage reconciliation cascade is:
#       1. Exact match
#       2. Normalised match  (case, whitespace, underscores)
#       3. Synonym resolution via GNVerifier (authority = "gnverifier", httr2)
#       4. Fuzzy match
#     Outputs saved to outputs/phylogeny/:
#       proceed_tree.rds        — pruned, dated phylogenetic tree
#       proceed_A_matrix.rds    — correlation matrix (from pr_phylo_cor)
#       proceed_name_map.rds    — data.frame mapping original sp_ncbi → canonical name
#       unmatched_species.csv   — species that remain unresolved after all stages
#
#   Stage 2 — apply_phylo_name_map()
#     Applies the saved name-mapping table to a data frame, adding a new column
#     sp_ncbi_canonical that holds the canonical tree-tip label.  This is done
#     once to the effect-size dataset in chapter 05 BEFORE entering the model
#     loop, so that every model sees sp_ncbi values that align with rownames(A).
#
#   Stage 3 — prepare_phylo_and_data()
#     For each model, subsets / expands A to the species in dat_model, filters
#     out NA rows, and builds V.  Called inside the chapter 05 loop.
# ─────────────────────────────────────────────────────────────────────────────

# ── Rotl fallback: direct tree retrieval with pruned-OTT-ID filtering ─────────
#
# The Open Tree of Life API returns HTTP 400 when tol_induced_subtree() is
# called with an OTT ID that has been pruned from the current synthetic tree.
# pr_get_tree() does not pre-filter these.
#
# This helper calls rotl::is_in_tree() on every matched OTT ID before building
# the subtree.  is_in_tree() returns FALSE for pruned IDs, so they are removed
# before the tol_induced_subtree() call.
#
# Note: tol_induced_subtree() returns tip labels as "Genus_species_ott1234567".
# .clean_ott_tip_labels() strips the OTT suffix and converts underscores to
# spaces so that reconcile_tree() receives plain binomial names, matching what
# pr_get_tree() delivers when it succeeds.
.clean_ott_tip_labels <- function(tree) {
  tree$tip.label <- sub("_ott[0-9]+$", "", tree$tip.label)
  tree$tip.label <- gsub("_", " ", tree$tip.label)
  tree
}

.build_tree_rotl_direct <- function(species_vec) {
  # context_name is intentionally omitted so tnrs uses "All life" — the PROCEED
  # database spans many taxa (fish, birds, plants, invertebrates, etc.), and
  # restricting to "Animals" would silently miss non-animal species.
  taxa <- tryCatch(
    rotl::tnrs_match_names(species_vec),
    error = function(e) {
      message(".build_tree_rotl_direct — tnrs_match_names failed: ",
              conditionMessage(e))
      NULL
    }
  )
  if (is.null(taxa) || nrow(taxa) == 0) return(NULL)

  taxa_ok <- taxa[!is.na(taxa$ott_id), ]
  if (nrow(taxa_ok) == 0) {
    message(".build_tree_rotl_direct — no OTT IDs matched.")
    return(NULL)
  }

  # is_in_tree() returns a logical vector — FALSE for pruned_ott_id entries
  in_tree <- tryCatch(
    rotl::is_in_tree(taxa_ok$ott_id),
    error = function(e) {
      message(".build_tree_rotl_direct — is_in_tree() failed: ",
              conditionMessage(e),
              "\nAssuming all IDs are in tree.")
      rep(TRUE, nrow(taxa_ok))
    }
  )

  n_pruned     <- sum(!in_tree)
  taxa_in_tree <- taxa_ok[in_tree, ]

  if (n_pruned > 0) {
    message("  is_in_tree() removed ", n_pruned,
            " pruned OTT ID(s): ",
            paste(taxa_ok$unique_name[!in_tree], collapse = ", "))
  }

  if (nrow(taxa_in_tree) < 3) {
    message("  Too few species remain in the synthetic tree (",
            nrow(taxa_in_tree), ") — cannot build subtree.")
    return(NULL)
  }

  tree_out <- tryCatch(
    rotl::tol_induced_subtree(ott_ids = taxa_in_tree$ott_id),
    error = function(e) {
      message(".build_tree_rotl_direct — tol_induced_subtree failed: ",
              conditionMessage(e))
      NULL
    }
  )

  # Clean OTT-format tip labels ("Genus_species_ott1234567" → "Genus species")
  # so that reconcile_tree() receives plain binomial names.
  if (!is.null(tree_out)) .clean_ott_tip_labels(tree_out) else NULL
}


# ── Stage 1: build or read ────────────────────────────────────────────────────

build_or_read_phylogeny <- function(
    species_vec,
    tree_path      = here::here("Rdata", "phylogeny", "proceed_tree.rds"),
    A_path         = here::here("Rdata", "phylogeny", "proceed_A_matrix.rds"),
    name_map_path  = here::here("Rdata", "phylogeny", "proceed_name_map.rds"),
    unmatched_path = here::here("Rdata", "phylogeny", "unmatched_species.csv"),
    recompute      = FALSE) {

  if (file.exists(A_path) && file.exists(name_map_path) && !recompute) {
    message("Reading cached phylogenetic matrix: ", A_path)
    return(list(
      A        = readRDS(A_path),
      name_map = readRDS(name_map_path)
    ))
  }

  # Read-only book render: never rebuild the phylogeny (network calls) on the fly.
  if (isTRUE(getOption("pace.read_only", FALSE)) && !recompute) {
    stop("[read-only render] Phylogeny cache not found (",
         A_path, " / ", name_map_path,
         ").\nThe book will not rebuild the tree. Sync the .rds files, or set recompute = TRUE deliberately.",
         call. = FALSE)
  }

  message("Building phylogenetic scaffold via prepR4pcm ...")

  # ── 1. Collect unique, non-empty species names ──────────────────────────────
  sp_raw <- unique(as.character(species_vec))
  sp_raw <- sp_raw[!is.na(sp_raw) & nzchar(sp_raw)]

  if (length(sp_raw) < 3) {
    message("Too few species (", length(sp_raw), ") to build a tree.")
    return(NULL)
  }

  # ── 2. Normalise names (handles underscores, case, extra whitespace) ─────────
  sp_norm <- tryCatch(
    prepR4pcm::pr_normalize_names(sp_raw),
    error = function(e) {
      message("pr_normalize_names failed: ", conditionMessage(e),
              "\nFalling back to raw names.")
      sp_raw
    }
  )

  # ── 3. Fetch tree via rotl ────────────────────────────────────────────────────
  # Try pr_get_tree() first.  If the OTL API returns a 400 because one or more
  # OTT IDs are flagged as "pruned_ott_id" in the synthetic tree, catch that
  # error and fall back to a direct rotl call that filters with is_in_tree()
  # before calling tol_induced_subtree().
  message("Fetching tree from Open Tree of Life ...")
  tree_raw <- tryCatch(
    prepR4pcm::pr_get_tree(sp_norm, source = "rotl"),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("pruned_ott_id", msg, fixed = TRUE) ||
          grepl("HTTP failure: 400", msg, fixed = TRUE)) {
        message(
          "pr_get_tree hit pruned OTT IDs (HTTP 400).\n",
          "Retrying with direct rotl + is_in_tree() filter ..."
        )
        .build_tree_rotl_direct(sp_norm)
      } else {
        message("pr_get_tree failed: ", msg)
        NULL
      }
    }
  )

  if (is.null(tree_raw)) {
    message("Tree retrieval failed. Returning NULL.")
    return(NULL)
  }

  # pr_get_tree() returns a pr_tree_result list; reconcile_tree() expects a
  # plain ape::phylo.  Extract $tree when needed.
  # The fallback (.build_tree_rotl_direct) already returns a plain phylo.
  tree_phylo <- if (inherits(tree_raw, "pr_tree_result")) tree_raw$tree else tree_raw

  if (is.null(tree_phylo) || !inherits(tree_phylo, "phylo")) {
    message("Tree object is not a valid ape::phylo — cannot proceed.")
    return(NULL)
  }

  # ── 4. Reconcile species names against tree tip labels (4-stage cascade) ─────
  # Synonym resolution uses the GNVerifier web service (authority = "gnverifier",
  # via httr2) rather than the local taxadb package — gnverifier resolves far
  # more names (needs no local taxonomic database) and requires network access.
  message("Reconciling species names (exact → normalised → synonym → fuzzy) ...")
  dat_temp <- data.frame(sp_ncbi = sp_norm, stringsAsFactors = FALSE)

  recon <- tryCatch(
    prepR4pcm::reconcile_tree(
      dat_temp,
      tree_phylo,
      x_species = "sp_ncbi",   # correct param (not species_col)
      fuzzy     = TRUE,
      authority = "gnverifier"
    ),
    error = function(e) {
      message("reconcile_tree failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(recon)) return(NULL)

  # ── 5. Apply reconciliation: prune tree to resolved species ─────────────────
  # drop_unresolved = TRUE prunes unmatched tree tips and removes unmatched data
  # rows from result$data.  The pruned tree has tip labels renamed to data-side
  # names (sp_norm values), which become the rownames of A after pr_phylo_cor().
  result <- tryCatch(
    prepR4pcm::reconcile_apply(recon,
                                data            = dat_temp,
                                tree            = tree_phylo,
                                species_col     = "sp_ncbi",
                                drop_unresolved = TRUE),
    error = function(e) {
      message("reconcile_apply failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(result)) return(NULL)

  tree_final <- result$tree

  # ── 6. Extract name-mapping table ───────────────────────────────────────────
  # reconcile_mapping() columns:
  #   name_x        — data-side name (sp_norm value we supplied)
  #   name_y        — original tree-tip label
  #   name_resolved — canonical name from taxonomic authority (may be NA)
  #   in_x / in_y   — logical presence flags
  #   match_type    — "exact" | "normalized" | "synonym" | "fuzzy" | "unresolved"
  #
  # After reconcile_apply(drop_unresolved=TRUE), tree tips are renamed to name_x
  # values, so rownames(A) == name_x for matched species.
  # We build: sp_raw → sp_ncbi_canonical (= name_x for resolved, NA otherwise).
  name_map_raw <- tryCatch(
    prepR4pcm::reconcile_mapping(recon),
    error = function(e) {
      message("reconcile_mapping failed: ", conditionMessage(e))
      NULL
    }
  )

  name_map <- if (!is.null(name_map_raw) &&
                   all(c("name_x", "in_x", "in_y") %in% names(name_map_raw))) {

    # Filter to data-side rows only (in_x = TRUE).
    # Tree-only rows (in_x = FALSE) have name_x = NA and would poison the lookup.
    data_rows <- name_map_raw[name_map_raw$in_x, ]

    # Build a lookup: sp_norm → canonical A-tip label.
    # Resolved species (in_x & in_y) get their sp_norm value as canonical (= tip label).
    # Unresolved species get NA → prepare_phylo_and_data() will assign star correlation.
    canonical_lookup <- stats::setNames(
      ifelse(data_rows$in_y, data_rows$name_x, NA_character_),
      data_rows$name_x
    )

    # Print reconciliation summary by match_type
    if ("match_type" %in% names(data_rows)) {
      resolved_rows <- data_rows[data_rows$in_y, ]
      mt <- table(resolved_rows$match_type)
      message("Name reconciliation: ",
              nrow(resolved_rows), "/", nrow(data_rows),
              " data species resolved.")
      for (m in names(mt)) message("  ", m, ": ", mt[[m]])
    }

    # Final table: one row per unique sp_raw value
    data.frame(
      sp_ncbi_original  = sp_raw,
      sp_ncbi_canonical = canonical_lookup[sp_norm],
      stringsAsFactors  = FALSE
    )

  } else {
    # Fallback: identity map from tree tip labels (no synonym resolution info)
    message("reconcile_mapping() result lacks expected columns; ",
            "using tree tip labels as canonical names.")
    data.frame(
      sp_ncbi_original  = tree_final$tip.label,
      sp_ncbi_canonical = tree_final$tip.label,
      stringsAsFactors  = FALSE
    )
  }

  # Save unmatched species for audit
  unmatched <- name_map[is.na(name_map$sp_ncbi_canonical) |
                           !nzchar(name_map$sp_ncbi_canonical), ]
  if (nrow(unmatched) > 0) {
    readr::write_csv(unmatched, unmatched_path)
    message(nrow(unmatched), " species unresolved — saved to ", unmatched_path)
  }

  # ── 7. Ensure branch lengths ─────────────────────────────────────────────────
  tree_final <- ape::multi2di(tree_final)
  if (is.null(tree_final$edge.length) || all(is.na(tree_final$edge.length))) {
    message("No branch lengths — assigning Grafen lengths.")
    tree_final <- ape::compute.brlen(tree_final, method = "Grafen", power = 1)
  }

  # ── 8. Build phylogenetic correlation matrix ─────────────────────────────────
  A <- tryCatch(
    prepR4pcm::pr_phylo_cor(tree_final),
    error = function(e) {
      message("pr_phylo_cor failed (", conditionMessage(e), "); using ape::vcv.phylo.")
      ape::vcv.phylo(tree_final, corr = TRUE)
    }
  )

  # ── 9. Save outputs ───────────────────────────────────────────────────────────
  saveRDS(tree_final, tree_path)
  saveRDS(A,          A_path)
  saveRDS(name_map,   name_map_path)

  message("Phylogenetic matrix saved: ", nrow(A), " × ", ncol(A),
          " (", sum(!is.na(name_map$sp_ncbi_canonical)), "/",
          length(sp_raw), " species resolved)")

  list(A = A, name_map = name_map)
}


# ── Stage 2: apply name map to dataset ───────────────────────────────────────
#
# Adds sp_ncbi_canonical column: the canonical tree-tip name for each row.
# Rows whose sp_ncbi did not resolve to a tree tip get NA in sp_ncbi_canonical
# (these are handled transparently by prepare_phylo_and_data below).
#
# Call this ONCE on dat_es in chapter 05, before entering the model loop.
apply_phylo_name_map <- function(dat,
                                  name_map,
                                  original_col  = "sp_ncbi",
                                  canonical_col = "sp_ncbi_canonical") {
  if (is.null(name_map)) {
    dat[[canonical_col]] <- dat[[original_col]]
    return(dat)
  }

  lookup <- stats::setNames(name_map$sp_ncbi_canonical,
                             name_map$sp_ncbi_original)

  dat[[canonical_col]] <- lookup[as.character(dat[[original_col]])]
  dat
}


# ── Stage 3: align data to A for a single model ──────────────────────────────
#
# Two sources of row loss, now clearly separated:
#   1. sp_ncbi_canonical is NA  — species had no resolved tree-tip name.
#      (Includes original NAs AND species the cascade could not match.)
#   2. sp_ncbi is NA / blank    — no species name at all in the original data.
#
# fill_unmatched = TRUE  (default): rows that resolved to a canonical name are
#   kept; A is expanded for any that ended up outside the tree using zero
#   off-diagonal (star phylogeny).  This is the primary-model behaviour.
#
# fill_unmatched = FALSE (strict): only rows whose canonical name appears as a
#   rowname in A_full are kept.  Use for the no-phylogeny sensitivity.
expand_A_with_unmatched <- function(A_full, all_species) {
  all_species <- all_species[!is.na(all_species) & nzchar(all_species)]
  already_in  <- all_species[all_species %in% rownames(A_full)]
  missing     <- setdiff(all_species, rownames(A_full))

  if (length(missing) == 0) {
    return(A_full[already_in, already_in, drop = FALSE])
  }

  n_all   <- length(all_species)
  A_exp   <- matrix(0, nrow = n_all, ncol = n_all,
                    dimnames = list(all_species, all_species))
  diag(A_exp) <- 1

  if (length(already_in) >= 2) {
    A_exp[already_in, already_in] <- A_full[already_in, already_in]
  }

  A_exp
}

prepare_phylo_and_data <- function(dat_model, A_full,
                                   fill_unmatched = TRUE,
                                   min_rows       = 10,
                                   min_species    = 2,
                                   label          = "") {
  pfx <- if (nzchar(label)) paste0("[", label, "] ") else ""

  has_phylo      <- !is.null(A_full)
  A_mod          <- NULL
  n_excl_na      <- 0L
  n_excl_noname  <- 0L
  n_unmatched_sp <- 0L

  # Decide which species column to use for matching A:
  # prefer the canonical column if it was created by apply_phylo_name_map()
  sp_col <- if ("sp_ncbi_canonical" %in% names(dat_model)) {
    "sp_ncbi_canonical"
  } else {
    "sp_ncbi"
  }

  # ── Step 1: Remove rows with no usable species name ─────────────────────────
  n_before  <- nrow(dat_model)
  dat_model <- dat_model |>
    dplyr::filter(!is.na(.data[[sp_col]]) & nzchar(as.character(.data[[sp_col]])))
  n_excl_na <- n_before - nrow(dat_model)

  if (n_excl_na > 0) {
    message(pfx, n_excl_na,
            " rows excluded: ", sp_col, " is NA or unresolved ",
            "(no canonical species name available).")
  }

  if (!has_phylo || nrow(dat_model) < min_rows) {
    dat_model <- dat_model |>
      dplyr::mutate(es_id_model = factor(seq_len(dplyr::n())))
    V <- diag(dat_model$vi_lnM_safe)
    rownames(V) <- colnames(V) <- as.character(dat_model$es_id_model)
    return(list(dat_model = dat_model, A_mod = NULL, V = V,
                has_phylo = FALSE,
                n_excl_na = n_excl_na, n_excl_noname = 0L,
                n_unmatched_sp = 0L))
  }

  # ── Step 2: Build (possibly expanded) A ─────────────────────────────────────
  all_species <- unique(as.character(dat_model[[sp_col]]))
  sp_in_tree  <- all_species[all_species %in% rownames(A_full)]
  sp_missing  <- setdiff(all_species, rownames(A_full))

  if (fill_unmatched) {
    if (length(sp_missing) > 0) {
      message(pfx, length(sp_missing),
              " species absent from A — assigned star (zero) correlation. ",
              "Examples: ",
              paste(head(sp_missing, 5), collapse = ", "),
              if (length(sp_missing) > 5) ", ..." else "")
    }
    n_unmatched_sp <- length(sp_missing)

    if (length(sp_in_tree) >= min_species) {
      A_mod     <- expand_A_with_unmatched(A_full, all_species)
      has_phylo <- TRUE
    } else {
      message(pfx, "Too few tree-matched species (", length(sp_in_tree),
              "); fitting without phylogenetic term.")
      has_phylo <- FALSE
    }

    dat_model <- dat_model |>
      dplyr::mutate(
        sp_ncbi     = droplevels(factor(as.character(.data[[sp_col]]),
                                        levels = all_species)),
        es_id_model = factor(seq_len(dplyr::n()))
      )

  } else {
    # Strict mode: exclude species not in A_full
    if (length(sp_in_tree) >= min_species) {
      A_mod <- A_full[sp_in_tree, sp_in_tree, drop = FALSE]

      dat_strict    <- dat_model |>
        dplyr::filter(as.character(.data[[sp_col]]) %in% sp_in_tree) |>
        dplyr::mutate(
          sp_ncbi     = droplevels(factor(as.character(.data[[sp_col]]),
                                          levels = sp_in_tree)),
          es_id_model = factor(seq_len(dplyr::n()))
        )
      n_excl_noname <- nrow(dat_model) - nrow(dat_strict)

      if (n_excl_noname > 0) {
        message(pfx, n_excl_noname,
                " rows excluded: species not in A (strict mode).")
      }

      if (nrow(dat_strict) >= min_rows) {
        dat_model <- dat_strict
        has_phylo <- TRUE
      } else {
        message(pfx, "Too few rows (", nrow(dat_strict),
                ") after strict filter; dropping phylogenetic term.")
        has_phylo <- FALSE
        A_mod     <- NULL
        dat_model <- dat_model |>
          dplyr::mutate(
            sp_ncbi     = droplevels(factor(as.character(.data[[sp_col]]))),
            es_id_model = factor(seq_len(dplyr::n()))
          )
      }
    } else {
      has_phylo <- FALSE
      message(pfx, "Only ", length(sp_in_tree),
              " species in A (strict mode); dropping phylogenetic term.")
      dat_model <- dat_model |>
        dplyr::mutate(es_id_model = factor(seq_len(dplyr::n())))
    }
  }

  # ── Step 3: Build V ──────────────────────────────────────────────────────────
  V <- diag(dat_model$vi_lnM_safe)
  rownames(V) <- colnames(V) <- as.character(dat_model$es_id_model)

  list(
    dat_model      = dat_model,
    A_mod          = A_mod,
    V              = V,
    has_phylo      = has_phylo,
    n_excl_na      = n_excl_na,
    n_excl_noname  = n_excl_noname,
    n_unmatched_sp = n_unmatched_sp
  )
}
