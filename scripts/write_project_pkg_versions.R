# scripts/write_project_pkg_versions.R
# Generate a Markdown table of the project packages (from install_pkgs.R)
# and insert it into README.md between markers: <!-- PKG-VERSIONS-START --> and <!-- PKG-VERSIONS-END -->

trim <- function(x) gsub("^\s+|\s+$", "", x)

read_required_pkgs <- function(path = "install_pkgs.R") {
  if (!file.exists(path)) stop("install_pkgs.R not found at ", path)
  txt <- readLines(path, warn = FALSE)
  start <- grep("required_pkgs\\s*<-\\s*c\\(", txt)
  if (length(start) == 0) stop("required_pkgs <- c(...) not found in install_pkgs.R")
  start <- start[1]
  rem <- txt[start:length(txt)]
  # find the line with the closing parenthesis for the c(...)
  rel_end <- which(grepl("\\)", rem))[1]
  if (is.na(rel_end)) stop("Could not find end of required_pkgs vector")
  frag <- paste(rem[1:rel_end], collapse = "\n")
  inner <- sub(".*c\\(", "", frag)
  inner <- sub("\\).*", "", inner)
  pkgs <- unlist(strsplit(inner, ","))
  pkgs <- trim(gsub("[\"']", "", pkgs))
  pkgs[pkgs != ""]
}

read_pinned_versions <- function(path = "install_pkgs.R") {
  txt <- readLines(path, warn = FALSE)
  pats <- gregexpr("remotes::install_version\\s*\\(\\s*\"([^\"]+)\"\\s*,\\s*version\\s*=\\s*\"([^\"]+)\"", txt, perl = TRUE)
  pins <- list()
  for (i in seq_along(pats)) {
    m <- pats[[i]]
    if (m[1] != -1) {
      matches <- regmatches(txt[i], m)
      for (mm in matches) {
        # mm like remotes::install_version("RcppArmadillo", version = "0.12.8.4.0")
        pkg <- sub(".*remotes::install_version\\(\\s*\"([^\"]+)\".*", "\\1", mm)
        ver <- sub(".*version\\s*=\\s*\"([^\"]+)\".*", "\\1", mm)
        pins[[pkg]] <- ver
      }
    }
  }
  pins
}

read_local_tarballs <- function(path = "install_pkgs.R") {
  txt <- readLines(path, warn = FALSE)
  # find install.packages(... "external_codes/*.tar.gz" ...)
  lines <- grep("external_codes/.*\\.tar\\.gz", txt, value = TRUE)
  tars <- character(0)
  for (ln in lines) {
    m <- regmatches(ln, regexpr("external_codes/[^\)\"]+\\.tar\\.gz", ln))
    if (length(m) && nchar(m)) tars <- c(tars, m)
  }
  # return named vector of file-> guessed pkg/version
  res <- list()
  for (t in unique(tars)) {
    f <- t
    bn <- basename(f)
    bn_noext <- sub("\\.tar\\.gz$", "", bn)
    # try split at last _ to separate pkg and version
    last_underscore <- max(regexpr("_", bn_noext), 0)
    if (last_underscore > 0) {
      name <- substr(bn_noext, 1, last_underscore-1)
      ver <- substr(bn_noext, last_underscore+1, nchar(bn_noext))
    } else {
      name <- bn_noext
      ver <- ""
    }
    res[[name]] <- list(file = f, guessed_version = ver)
  }
  res
}

make_table <- function(pkgs, pins, tars) {
  inst <- as.data.frame(installed.packages()[, c("Package", "Version")], stringsAsFactors = FALSE)
  rownames(inst) <- inst$Package
  rows <- list()
  for (p in pkgs) {
    installed_v <- if (p %in% rownames(inst)) inst[p, "Version"] else "NOT INSTALLED"
    pinned_v <- if (!is.null(pins[[p]])) pins[[p]] else ""
    rows[[p]] <- c(p, pinned_v, installed_v)
  }
  # include pinned local tarballs if their name matches not in pkgs
  for (name in names(tars)) {
    if (!(name %in% pkgs)) {
      installed_v <- if (name %in% rownames(inst)) inst[name, "Version"] else "NOT INSTALLED"
      pinned_v <- tars[[name]]$guessed_version
      rows[[name]] <- c(name, pinned_v, installed_v)
    }
  }
  # build md
  md <- c()
  md <- c(md, "| Package | Pinned / Required | Installed |")
  md <- c(md, "|---|---|---|")
  for (r in rows) {
    md <- c(md, sprintf("| %s | %s | %s |", r[1], ifelse(nzchar(r[2]), r[2], "-"), r[3]))
  }
  md
}

insert_into_readme <- function(md_lines, readme_path = "README.md") {
  if (!file.exists(readme_path)) stop("README.md not found")
  rd <- readLines(readme_path, warn = FALSE)
  start <- grep("<!-- PKG-VERSIONS-START -->", rd)
  end <- grep("<!-- PKG-VERSIONS-END -->", rd)
  block <- c("<!-- PKG-VERSIONS-START -->", "", md_lines, "", "<!-- PKG-VERSIONS-END -->")
  if (length(start) == 1 && length(end) == 1 && start < end) {
    rd2 <- c(rd[1:(start-1)], block, rd[(end+1):length(rd)])
  } else {
    # append to end
    rd2 <- c(rd, "", block)
  }
  writeLines(rd2, readme_path)
  message("README.md updated with package versions table between markers.")
}

main <- function() {
  pkgs <- tryCatch(read_required_pkgs(), error = function(e) {
    message("Warning: could not extract required_pkgs: ", e$message)
    character(0)
  })
  pins <- tryCatch(read_pinned_versions(), error = function(e) list())
  tars <- tryCatch(read_local_tarballs(), error = function(e) list())
  md_lines <- make_table(pkgs, pins, tars)
  # add header with R version and timestamp
  header <- sprintf("**R version:** %s  ", R.version.string)
  header2 <- sprintf("**Generated:** %s", format(Sys.time(), tz = Sys.timezone()))
  insert_into_readme(c(header, header2, "", md_lines))
}

if (sys.nframe() == 0) main()
