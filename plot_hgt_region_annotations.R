setwd("/Research/postdoc/p_apium/hgt")
##read in HGT file##
hgt <- read.table("hgt_region_genes.txt", sep = "\t", header = TRUE)

##read in annotations file
annotations <- read.table("hgt_annotations.txt", sep = "\t", header = TRUE)

##variables for high and low xlim values##
xlows <- c(393079, 1227979, 1464686, 1583893, 1847636)
xhighs <- c(410932, 1238433, 1470341, 1589872, 1858483)

elongation <- 197 ##length of the arrow on the gene boxes
par(mar = c(1.1, 1.1, 1.1, 1.1))
i <- 5
##for (i in unique(hgt$hgt_region)) { ##start looping through HGT dataframe
  ##variable for use in adding domain labels
  positive_strand_counter <- 0

  ##make pdf##
##  pdf(paste("hgt_region_", i, "_detail.pdf", sep = "")) #, width = 10.5, height = 8)
  svg(paste("hgt_region_", i, "_detail.svg", sep = "")) #, width = 10.5, height = 8)
  
  ##plot window for HGT region ##
  plot(-1000, -1000, xlim = c(xlows[i], xhighs[i]), ylim = c(0,1), xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
  
  ##draw representation of genome region##
  abline(h = 0.5)

  for (k in 1:length(hgt$start[hgt$hgt_region == i])) { ##loop through each gene in the HGT region
    gene <- hgt$ortho_name[hgt$hgt_region == i][k] ##gene name for reference later
    gene_start <- hgt$start[hgt$ortho_name == gene] ##gene start
    gene_end <- hgt$end[hgt$ortho_name == gene] ##gene end
    
    ##add gene representations##
    if (hgt$strand[hgt$ortho_name == gene] == "neg") { ##if gene is on negative strand, draw underneath genome line with arrow pointing to the left
      polygon(c(gene_start + elongation, gene_start, gene_start + elongation, gene_end, gene_end), c(0.4, 0.45, 0.5, 0.5, 0.4), col = rgb(0.5, 0.5, 0.5, 0.25), border = "black")
    } else { ##if gene is on positive strand, draw above genome line with arrow pointing to the right
      polygon(c(gene_start, gene_start, gene_end - elongation, gene_end, gene_end - elongation), c(0.5, 0.6, 0.6, 0.55, 0.5), col = rgb(0.5, 0.5, 0.5, 0.25), border = "black")
    }
    
    ##add domains to gene representations##
    for (j in 1:length(annotations$ortho_name[annotations$ortho_name == gene])) { ##some genes have more than one domain annotation, loop through each annotation for the gene
      domain_start <- annotations$start[annotations$ortho_name == gene & annotations$annotation_num == j] + gene_start
      domain_end <- annotations$end[annotations$ortho_name == gene & annotations$annotation_num == j] + gene_start
      lab <- annotations$annotation[annotations$ortho_name == gene & annotations$annotation_num == j]
      if (hgt$strand[hgt$ortho_name == gene] == "neg") { ##if gene is on negative strand, draw underneath genome line
        polygon(c(domain_start, domain_start, domain_end, domain_end), c(0.4, 0.5, 0.5, 0.4), col = rgb(0.78125, 0.3984375, 1, 0.5), border = NA)
        if (gene == "paap|1085") { ##make adjustment so the domain name fits on the plot
          text(((domain_start + domain_end) / 2) + 500, 0.35, labels = lab, cex = 0.7)
        } else {
          text((domain_start + domain_end) / 2, 0.35, labels = lab, cex = 0.7)
        }
      } else { ##if gene is on positive strand, draw above genome line
        if (positive_strand_counter %% 2 == 0) { ##alternate writing names above and below domains, and alternate colors
          y <- 0.65
          red <- 0.7843137
          green <- 0.4
          blue = 1
        } else {
          y <- 0.45
          red <- 1
          green <- 0.4
          blue <- 0.4 
        }
        polygon(c(domain_start, domain_start, domain_end, domain_end), c(0.5, 0.6, 0.6, 0.5), col = rgb(red, green, blue, 0.5), border = NA)
        if (gene == "paap|1283" & j == 1) { ##make adjustment so the domain name fits on the plot
          text(((domain_start + domain_end) / 2) + 300 , y + 0.025, labels = lab, cex = 0.7)
        } else {
          text((domain_start + domain_end) / 2, y, labels = lab, cex = 0.7)          
        }
        positive_strand_counter <- positive_strand_counter + 1
      }
    }
  } 
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  dev.off()
  ##}
