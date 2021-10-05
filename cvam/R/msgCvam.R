.msgCvam <- function( msg.codes, msg.len.actual ){
   ###########################################
   # private function 
   # converts matrix of integer message codes
   # to character text;
   # relies on internal package data object
   # icodesCvam
   ###########################################
   if( msg.len.actual == 0 ){
      msg.lines <- NULL
   }
   else{
      msg.lines <- character( msg.len.actual )
      for( i in 1:msg.len.actual ){
         if( msg.codes[i,1L] == 1L ){
            code <- msg.codes[i,4L]
            if( ( code > 0 ) & ( code <=  length( icodesCvam$comments ) ) ){
               msg.lines[i] <- icodesCvam$comments[ code ]
            }
            else{
               msg.lines[i] <- "???"
            }
         }
         else if( msg.codes[i,1L] == 2L ){
            code <- msg.codes[i,3L]
            if( ( code > 0 ) & ( code <= length(icodesCvam$subnames ) ) ){
               subname <- icodesCvam$subnames[ code ]
            }
            else{
               subname <- "???"
            }
            code <- msg.codes[i,2L]
            if( ( code > 0 ) & ( code <= length(icodesCvam$modnames ) ) ){
               modname <- icodesCvam$modnames[ code ]
            }
            else{
               modname <- "???"
            }
            msg.lines[i] <- paste( "OCCURRED IN:", subname, "in MOD",
               modname, sep=" " )
         }
         else if( msg.codes[i,1L] == 3L ){
            msg.lines[i] <- paste( "Observation", 
   	    format( msg.codes[i,5L] ),
               sep=" " )
         }
         else if( msg.codes[i,1L] == 4L ){
            msg.lines[i] <- paste( "Variable", 
   	    format( msg.codes[i,6L] ),
               sep=" " )
         }
         else if( msg.codes[i,1L] == 5L ){
            msg.lines[i] <- paste( "Iteration", 
   	    format( msg.codes[i,7L] ),
               sep=" " )
         }
         else if( msg.codes[i,1L] == 6L ){
            msg.lines[i] <- paste( "Iteration ", 
   	    format( msg.codes[i,7L] ),
               ", Cycle ",
   	    format( msg.codes[i,8L] ),
               sep="" )
         }
         else if( msg.codes[i,1L] == 7L ){
            msg.lines[i] <- paste( "Group", 
   	    format( msg.codes[i,9L] ),
               sep=" " )
         }
         else if( msg.codes[i,1L] == 8L ){
            msg.lines[i] <- paste( "Group", 
   	    format( msg.codes[i,9L] ),
               ", Term ",
   	    format( msg.codes[i,10L] ),
               sep="" )
         }
         else if( msg.codes[i,1L] == 9L ){
            msg.lines[i] <- paste( "Factor", 
   	    format( msg.codes[i,11L] ),
               sep=" " )
         }
         else if( msg.codes[i,1L] == 10L ){
            msg.lines[i] <- paste( "Factor ", 
   	    format( msg.codes[i,11L] ),
               ", Level ",
   	    format( msg.codes[i,12L] ),
               sep="" )
         }
         else if( msg.codes[i,1L] == 11L ){
            msg.lines[i] <- paste( "Submodel", 
   	    format( msg.codes[i,13L] ),
               sep=" " )
         }
         else if( msg.codes[i,1L] == 12L ){
            msg.lines[i] <- paste( "Estimate", 
   	    format( msg.codes[i,14L] ),
               sep=" " )
         }
         else if( msg.codes[i,1L] == 13L ){
            msg.lines[i] <- paste( "Predict", 
   	    format( msg.codes[i,15L] ),
               sep=" " )
         }
         else if( msg.codes[i,1L] == 14L ){
            msg.lines[i] <- paste( "Impute", 
   	    format( msg.codes[i,16L] ),
               sep=" " )
         }
         else if( msg.codes[i,1L] == 15L ){
            msg.lines[i] <- paste( "Cell", 
   	    format( msg.codes[i,17L] ),
               sep=" " )
         }
         else{
            msg.lines[i] <- "???"
         }
      }
   }
   return( msg.lines ) }
