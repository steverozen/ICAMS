source("~/.Rprofile")

.workspace_home = getwd()

mysnapreview = function() {
  testthat::snapshot_review(
    path = file.path(.workspace_home, "tests/testtthat/")
  )
}

message(".Rproflie in ", getwd(), " was sourced")
