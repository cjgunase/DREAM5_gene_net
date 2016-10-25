get.weight.matrix <- function(expr.matrix, K="sqrt", nb.trees=10, input.idx=NULL, importance.measure="IncNodePurity", seed=NULL, trace=TRUE, ...) {
  # set random number generator seed if seed is given
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # to be nice, report when parameter importance.measure is not correctly spelled
  if (importance.measure != "IncNodePurity" && importance.measure != "%IncMSE") {
    stop("Parameter importance.measure must be \"IncNodePurity\" or \"%IncMSE\"")
  }
  # transpose expression matrix to (samples x genes)
  expr.matrix <- t(expr.matrix)
  # normalize expression matrix
  expr.matrix <- apply(expr.matrix, 2, function(x) { (x - mean(x)) / sd(x) } )
  # setup weight matrix
  num.samples <- dim(expr.matrix)[1]
  num.genes <- dim(expr.matrix)[2]
  gene.names <- colnames(expr.matrix)
  weight.matrix <- matrix(0.0, nrow=num.genes, ncol=num.genes)
  rownames(weight.matrix) <- gene.names
  colnames(weight.matrix) <- gene.names
  # get number of input genes, names of input genes
  if (is.null(input.idx)) {
    num.input.genes <- num.genes
    input.gene.names <- gene.names
  } else {
    num.input.genes <- length(input.idx)
    # input gene indices given as integers
    if (is.numeric(input.idx)) {
      input.gene.names <- gene.names[input.idx]
      # input gene indices given as names
    } else {
      input.gene.names <- input.idx
      # for security, abort if some input gene name is not in gene names
      missing.gene.names <- setdiff(input.gene.names, gene.names)
      if (length(missing.gene.names) != 0) {
        for (missing.gene.name in missing.gene.names) {
          cat(paste("Gene ", missing.gene.name,
                    " was not in the expression matrix\n", sep=""))
        }
        stop("Aborting computation")
      }
    }
  }
  # set mtry
  if (class(K) == "numeric") {
    mtry <- K
  } else if (K == "sqrt") {
    mtry <- round(sqrt(num.input.genes))
  } else if (K == "all") {
    mtry <- num.input.genes-1
  } else {
    stop("Parameter K must be \"sqrt\", or \"all\", or an integer")
  }
  if (trace) {
    cat(paste("Starting RF computations with ", nb.trees,
              " trees/target gene,\nand ", mtry,
              " candidate input genes/tree node\n",
              sep=""))
    flush.console()
  }
  # compute importances for every target gene
  for (target.gene.idx in seq(from=1, to=num.genes)) {
    if (trace) {
      cat(paste("Computing gene ", target.gene.idx, "/", num.genes, "\n", sep=""))
      flush.console()
    }
    target.gene.name <- gene.names[target.gene.idx]
    # remove target gene from input genes
    these.input.gene.names <- setdiff(input.gene.names, target.gene.name)
    x <- expr.matrix[,these.input.gene.names]
    y <- expr.matrix[,target.gene.name]
    rf <- randomForest(x, y, mtry=mtry, ntree=nb.trees, importance=TRUE, ...)
    im <- importance(rf)[,importance.measure]
    im.names <- names(im)
    weight.matrix[im.names, target.gene.name] <- im
  }
  return(weight.matrix / num.samples)
}       
# get.link.list: get sorted list of regulatory links (most likely link first)
#
# Parameters (required):
#    -- weight.matrix: weighted adjacency matrix as returned by get.weight.matrix
#
# Parameters (optional):
#    -- report.max: maximum number of links to report (default all links)
#
# Returns:
#    list of links in data frame. Each line of the data frame has format
#    regulatory_gene target_gene importance_score
#
get.link.list <- function(weight.matrix, report.max=NULL) {
  # set negative weights (for permutation of OOB importance) to 0.0 
  # weight.matrix[weight.matrix < 0.0] <- 0.0
  num.genes <- dim(weight.matrix)[1]
  genes <- colnames(weight.matrix)
  matrix.length <- length(weight.matrix)
  list.length <- num.genes * (num.genes - 1)
  if (!is.null(report.max) && report.max < list.length) {
    list.length <- report.max
  }
  # setup link list
  link.list <- data.frame(from.gene=rep("", list.length),
                          to.gene=rep("", list.length),
                          im=rep(0.0, list.length),
                          stringsAsFactors=FALSE)
  sorted.indices <- order(weight.matrix, decreasing=TRUE)
  # fill in link list
  index.number <- 1
  link.number <- 1
  while (index.number <= matrix.length && link.number <= list.length) {
    i <- sorted.indices[index.number]
    im <- weight.matrix[i]
    row.col <- lin.to.square(i, num.genes)
    row <- row.col[1]
    col <- row.col[2]
    # Only process weights off-diagonal
    if (row != col) {
      from.gene <- genes[row]
      to.gene <- genes[col]
      link.list[link.number,] <- list(from.gene, to.gene, im)
      link.number <- link.number + 1
    }
    index.number <-index.number + 1
  }
  return(link.list)
}

# load required packages
tryCatch( suppressMessages(library(randomForest)),
          error=function(e) { cat("Error: package randomForest must be installed\n");
            cat("Use install.packages(\"randomForest\")\n") })

# utility function to convert linear index to (row,col) index for matrix
lin.to.square <- function(i, nrow) {
  col <- ((i - 1) %/% nrow) + 1
  row <- ((i - 1) %% nrow) + 1
  return(c(row, col))
}

func <- function(i,topTFs){
  temp <- NULL
  yVal <- gene.exp[,i]
  yName <- names(gene.exp)[i]
  cv <-cv.spls(tf.exp,yVal,eta = seq(0.5,0.9,0.1),K = c(5:10), kappa=0.5, select="pls2", fit="simpls",scale.x=TRUE, scale.y=TRUE, plot.it=F);
  f <- spls(tf.exp,yVal,eta =cv$eta.opt,K = cv$K.opt);
  #f <- spls(tf.exp,yVal,eta =0.7,K = 10)
  coef.f <- coef(f)
  temp <- as.data.frame(coef.f)
  temp$FROM <- rownames(temp)
  temp <- temp[temp['V1'] > 0,]
  temp <- temp[ order(-temp[,1]),]
  temp['TO'] <- rep(yName, times = dim(temp)[1])
  temp <- temp[1:topTFs,]#top TFs for each Target Gene is seleted. if it has less than 20 NA fill fill
  temp <- na.omit(temp) #this will remove if NAs present
  if(length(temp$FROM) > 3){
    myvars <- names(tf.exp) %in% temp$FROM
    selectedTF <- tf.exp[myvars]
    selectedTF[yName] <- yVal
    selectedTF <- t(selectedTF)
    temp <- as.data.frame(get.link.list(get.weight.matrix(selectedTF,input.idx = as.character(temp$FROM)),report.max =topTFs))
  }else{
    temp <- NULL
  }
  
  
  return(temp)
}

