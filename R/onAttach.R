".onAttach" <- function(lib, pkg) {
  mylib <- dirname(system.file(package = pkg))
  title <- utils::packageDescription(pkg, lib.loc = mylib)$Title
  ver <- utils::packageDescription(pkg, lib.loc = mylib)$Version
  author <- utils::packageDescription(pkg, lib.loc = mylib)$Author
  packageStartupMessage(pkg, ": ", title, "\nVersion: ", ver, "\nAuthors: ", author, "\n")
}
