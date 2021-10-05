########################################################################
# functions for coarsened factors
#
# documented in
#    baseLevels.Rd
#    coarsened.Rd
#    dropCoarseLevels.Rd
#    is.naCoarsened
#    rep.coarsened
#
# todo (not high priority)
#    function for converting all coarseLevels to baseLevels, resulting is
#       a factor; maybe make this a non-default option for
#       dropCoarseLevels?
#    function for adding a coarse level to a coarsened factor
#    baseLevels<-  (changes the character name of one or more
#       baseLevels, throwing an error if the result would create
#       duplicates; if the value argument is of a different length,
#       then we need to drop some or add some
#       Would it be better to just reform the result by a call
#       to factor? Maybe.
#    coarseLevels<- (changes the character name of one or more 
#       coarseLevels
#    droplevels.coarsened: does not turn it into a factor, but 
#       keeps it as a coarsened factor
#         *   drop unused baseLevels
#         *   drop unused coarseLevels, except 
#         *   drop redundant coarseLevels? (no; if we want to do this,
#              we should explicitly recode them first,
#              making the redundant ones empty, then drop the unused levels.
#              Or maybe have this as an option, with a default of not
#              doing it.
#         *    drop levels that are being used, if they are in exclude,
#              setting the corresponding observations to NA
#         *    But we cannot drop NA
########################################################################

coarsened <- function(obj, levelsList=list(), warnIfCoarsened=TRUE ){
   if( is.coarsened(obj) ) {
      if( warnIfCoarsened ) {
         warning( gettext("Argument 'obj' is already a coarsened factor",
            domain = NA))
         if( !is.null(levelsList) || length(levelsList) > 0L )
            warning( gettext("Argument 'levelsList' ignored", domain = NA))
      }
      return(obj)
   }
   if( ! inherits(obj, "factor") ) stop( gettext(
      "Argument 'obj' is not a factor"), domain = NA )
   if( ! ( typeof(levelsList) == "list" ) ) stop( gettext(
      "Argument 'levelsList' is not a list"), domain = NA )
   if( any( is.na( levels(obj) ) ) ) stop( gettext(
      "Argument 'obj' has NA level"), domain = NA )
   # identify coarseLevels and baseLevels
   stopifnot( ! anyDuplicated( names(levelsList) ) )
   coarseLevels <- if (length(levelsList) == 0L) 
      character(0) else names(levelsList)
   baseLevels <- levels(obj)[ ! levels(obj) %in% coarseLevels ]
   stopifnot( length(baseLevels) > 1L )
   obj <- factor(obj, levels=c(baseLevels, coarseLevels) )
   baseLevelCodes <- match( baseLevels, levels(obj) )
   coarseLevelCodes <- match( coarseLevels, levels(obj) )
   obj <- addNA(obj)
   coarseLevels <- c( coarseLevels, NA )
   coarseLevelCodes <- c( coarseLevelCodes, nlevels(obj) )
   # create mapping matrix
   mapping <- matrix( 0L, length(coarseLevels), length(baseLevels), 
      dimnames=list(coarseLevels, baseLevels) )
   for( i in seq_along(levelsList) ){
      stopifnot( all( levelsList[[i]] %in% baseLevels ) )
      stopifnot( length(levelsList[[i]]) > 1L )
      mapping[i, baseLevels %in% levelsList[[i]] ] <- 1L
   }
   mapping[ length(coarseLevels), ] <- 1L
   if( any( duplicated(mapping) ) ) warning( gettext(
      "There are multiple 'coarseLevels' with the same mapping" ),
      domain = NA )
   obj <- structure(obj, class = c("coarsened", class(obj) ),
      mapping = mapping,
      baseLevels = baseLevels, coarseLevels = coarseLevels,
      nBaseLevels = length(baseLevels), nCoarseLevels = length(coarseLevels),
      baseLevelCodes = baseLevelCodes,
      coarseLevelCodes = coarseLevelCodes,
      latent = FALSE )
   if( is.null( attr(obj,"contrasts") ) ) {
      contr <- if( is.ordered(obj) ) contr.poly( length(baseLevels) ) else
         contr.sum( length(baseLevels) )
   }
   else {
      contr <- contrasts(obj)
   }
   contr <- rbind(contr, matrix(0, length(coarseLevels), NCOL(contr) ) )
   rownames(contr) <- c(baseLevels, coarseLevels)
   attr(obj,"contrasts") <- contr
   obj
}

is.coarsened <- function(x){
   inherits(x, "coarsened")
}

print.coarsened <- function(x, quote = FALSE, max.levels = NULL,
   width = getOption("width"), ...){
   stopifnot( is.coarsened(x) )
   ord <- is.ordered(x)
   if (length(x) == 0L)
      cat("coarsened(0)\n", sep="")
   else{
      xx <- character(length(x))
      xx[] <- as.character(x)
      keepAttrs <- setdiff( names(attributes(x)),
         c("levels", "class", "mapping", "baseLevels",
           "coarseLevels", "nBaseLevels", "nCoarseLevels", 
           "baseLevelCodes", "coarseLevelCodes", "contrasts", "latent") )
      attributes(xx)[keepAttrs] <- attributes(x)[keepAttrs]
      print(xx, quote = quote, ...)
   }
   max2 <- maxl <- if (is.null(max.levels)) 
      TRUE
   else max.levels
   if (maxl) {
      n1 <- length(blev <- encodeString(attr(x, "baseLevels"),
         quote = ifelse(quote, "\"", "")))
      colsep1 <- if (ord) " < " else " "
      T1 <- "  Base levels: "
      if (is.logical(maxl)) 
         maxl <- {
            width1 <- width - (nchar(T1, "w") + 3L + 1L + 
               3L)
            lenl <- cumsum(nchar(blev, "w") + nchar(colsep1, 
               "w"))
            if (n1 <= 1L || lenl[n1] <= width1) 
               n1
            else max(1L, which.max(lenl > width1) - 1L)
         }
      drop1 <- n1 > maxl
      cat(if (drop1) 
         paste(format(n1), ""), T1, paste(if (drop1) 
         c(blev[1L:max(1, maxl - 1)], "...", if (maxl > 1) blev[n1])
      else blev, collapse = colsep1), "\n", sep = "")
   }
   if (max2) {
      n2 <- length(clev <- encodeString(attr(x, "coarseLevels"),
         quote = ifelse(quote, "\"", "")))
      colsep2 <- " "
      T2 <- "Coarse levels: "
      if (is.logical(max2)) 
         max2 <- {
            width2 <- width - (nchar(T2, "w") + 3L + 1L + 
               3L)
            len2 <- cumsum(nchar(clev, "w") + nchar(colsep2, 
               "w"))
            if (n2 <= 1L || len2[n2] <= width2) 
               n2
            else max(1L, which.max(len2 > width2) - 1L)
         }
      drop2 <- n2 > max2
      cat(if (drop2) 
         paste(format(n2), ""), T2, paste(if (drop2) 
         c(clev[1L:max(1, max2 - 1)], "...", if (max2 > 1) clev[n2])
      else clev, collapse = colsep2), "\n", sep = "")
      cat( paste("Mapping", 
         if( drop1 ) " (some levels not shown):" else ":",
         "\n", sep=""), sep="")
      print( if (drop1) attr(x, "mapping")[,1L:max(1, maxl - 1)]
         else attr(x, "mapping"), ... )
   }
   if( attr(x,"latent" ) ) message( gettext(
      "This coarsened factor is latent." ), domain=NA )
   if (!isTRUE(val <- .valid.factor(x))) 
      warning(val)
   invisible(x)
}

baseLevels <- function(x) {
   stopifnot( is.coarsened(x) )
   attr(x, "baseLevels")
}

coarseLevels <- function(x) {
   stopifnot( is.coarsened(x) )
   attr(x, "coarseLevels")
}

nBaseLevels <- function(x) {
   stopifnot( is.coarsened(x) )
   attr(x, "nBaseLevels")
}

nCoarseLevels <- function(x) {
   stopifnot( is.coarsened(x) )
   attr(x, "nCoarseLevels")
}

mapping <- function(x) {
   stopifnot( is.coarsened(x) )
   attr(x, "mapping")
}

.packedMapping <- function(x) {
   stopifnot( is.coarsened(x) )
   pM <- nlevels(x)
   for( i in seq_along( baseLevels(x) ) ) {
      pM <- c(pM, 1L, i )
   }
   map <- attr(x, "mapping")
   for( i in seq_along( coarseLevels(x) ) ) {
      pM <- c(pM, sum(map[i,]), (1:NCOL(map))[map[i,]==1L] )
   }
   pM
}

dropCoarseLevels <- function(x) {
   # turns coarsened object into a factor, dropping the coarseLevels
   # and changing their corresponding observations to NA
   stopifnot( inherits(x, "factor") )
   if( ! is.coarsened(x) ) return(x)
   contr <- contrasts(x)[-attr(x,"coarseLevelCodes"),,drop=FALSE] 
   bL <- baseLevels(x)
   x <- droplevels.factor( x, exclude = coarseLevels(x) )
   if( all(is.na(x)) ) levels(x) <- bL
   structure( x, mapping = NULL,
      baseLevels = NULL, coarseLevels = NULL, 
      nBaseLevels = NULL, nCoarseLevels = NULL, 
      baseLevelCodes = NULL, coarseLevelCodes = NULL,
      latent = NULL,
      class = class(x)[ class(x)!="coarsened" ],
      contrasts = contr)
}

is.naCoarsened <- function(x) {
   stopifnot( inherits(x, "coarsened" ) )
   iCodeNA <- seq_along( levels(x) )[ is.na( levels(x) ) ]
   unclass(x) == iCodeNA
}

droplevels.coarsened <- function( x, ... ) {
   stopifnot( inherits(x, "coarsened" ) )
   stop( gettext( 
      "'droplevels' cannot be used with coarsened factors"),
      domain = NA )
}

relevel.coarsened <- function( x, ... ) {
   stopifnot( inherits(x, "coarsened" ) )
   stop( gettext( "'relevel' cannot be used with coarsened factors" ),
      domain = NA )
}

reorder.coarsened <- function( x, ... ) {
   stopifnot( inherits(x, "coarsened" ) )
   stop( gettext( "'reorder' cannot be used with coarsened factors" ),
      domain = NA )
}

rep.coarsened <- function(x, ...) {
   y <- NextMethod()
   structure(y,
      levels = attr(x, "levels"),
      class = attr(x, "class"),
      mapping = attr(x, "mapping"),
      baseLevels = attr(x, "baseLevels"),
      coarseLevels = attr(x, "coarseLevels"),
      nBaseLevels = attr(x, "nBaseLevels"),
      nCoarseLevels = attr(x, "nCoarseLevels"),
      baseLevelCodes = attr(x, "baseLevelCodes"),
      coarseLevelCodes = attr(x, "coarseLevelCodes"),
      latent = attr(x, "latent"),
      contrasts = attr(x, "contrasts") )
}

latentFactor <- function( n, levels = 2L ) {
   n <- as.integer(n)[1L]
   stopifnot( n >= 0L )
   if( is.character(levels) ) {
      stopifnot( length(levels) >= 2L )
   } else {
      levels <- as.integer(levels)[1L]
      stopifnot( levels >= 2L )
      levels <- as.character(1:levels)
   }
   obj <- factor( rep(NA,n), levels = levels )
   obj <- coarsened(obj)
   attr( obj, "latent" ) <- TRUE
   obj
}

is.latentFactor <- function( x ) {
   if( ! inherits(x, "coarsened") ) FALSE else attr(x, "latent")
}

`[.coarsened` <- function (x, ... ) {
    y <- NextMethod("[")
    attributes(y) <- attributes(x)
    y
}

`[[.coarsened` <- function (x, ... ) {
    y <- NextMethod("[[")
    attributes(y) <- attributes(x)
    y
}

`[<-.coarsened` <- function (x, ..., value) {
    y <- NextMethod("[<-")
    w <- is.na(y)
    if( any(w) ) y[w] <- NA
    y
}

`[[<-.coarsened` <- function (x, ..., value) {
    y <- NextMethod("[[<-")
    w <- is.na(y)
    if( any(w) ) y[w] <- NA
    y
}

