# Install required Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")

pkgs <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggraph")
for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
        cat("Installing", p, "...\n")
        BiocManager::install(p, update = FALSE, ask = FALSE)
    } else {
        cat(p, "already installed\n")
    }
}
cat("Done\n")
