
Sys.setenv(TERM = "dumb") # will be unset by tinytest
tinytest:::.onLoad()
expect_false(getOption("tt.pr.color"))



options(tt.pr.color=NULL)
Sys.setenv(NO_COLOR=1L)
tinytest:::.onLoad()
expect_false(getOption("tt.pr.color"))


# reset option set by .onLoad()
options(tt.pr.color=NULL)


