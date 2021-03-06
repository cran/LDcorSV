LD.Measures <- function(donnees, V = NA, S = NA, data = "G", supinfo = FALSE, na.presence = TRUE) {
    
    loc1 <- vector("character", 0)
    loc2 <- vector("character", 0)
    r2 <- numeric(0)
    M.r2 <- NA
    r2v <- numeric(0)
    M.r2v <- NA
    r2s <- numeric(0)
    M.r2s <- NA
    r2vs <- numeric(0)
    M.r2vs <- NA
    MAF.loc1 <- numeric(0)
    heterofreq.loc1 <- numeric(0)
    NAfreq.loc1 <- numeric(0)
    MAF.loc2 <- numeric(0)
    heterofreq.loc2 <- numeric(0)
    NAfreq.loc2 <- numeric(0)
    RES <- numeric(0)
    
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {abs(x - round(x)) < tol}
    
    if(is.null(rownames(donnees)) | is.null(colnames(donnees)))
        return(print("ERROR: genotype data must have row and column names"))
    
    mark <- colnames(donnees)
    
    if(all(!is.wholenumber(donnees[!is.na(donnees)])))
        return(print("ERROR: haplotype or genotype data must be a numeric matrix with integer values"))
    
    if(data == "G")
        if(max(donnees[!is.na(donnees)]) > 2 | min(donnees[!is.na(donnees)]) < 0)
            return(print("ERROR: genotype data must be a numeric matrix with integer values between 0 and 2"))
    
    if(data == "H")
        if(max(donnees[!is.na(donnees)]) > 1 | min(donnees[!is.na(donnees)]) < 0)
            return(print("ERROR: haplotype must be a numeric matrix with integer values between 0 and 1"))
    
    if(length(unique.default(is.na(S))) == 2)
        return(print("ERROR: missing values are not allowed in the structure matrix"))
    
    if(length(unique.default(is.na(V))) == 2)
        return(print("ERROR: missing values are not allowed in the genetic variance-covariance matrix"))
    
    S.yes <- all(is.na(S))
    V.yes <- all(is.na(V))
    
    if (!S.yes) {
        if (max(S) > 1 | min(S) < 0)
            return(print("ERROR: the structure matrix must be a numeric matrix with values between 0 and 1"))
        if (is.null(colnames(S)) | is.null(rownames(S)))
            return(print("ERROR: structure matrix  must have row and column names"))
    }
    
    if (!V.yes) {
        if (is.null(colnames(V)) | is.null(rownames(V)))
            return(print("ERROR: genotypic variance-covariance matrix  must have row and column names"))
    }
    
    if (!S.yes & !V.yes) {
        ID    <- rownames(donnees)
        ID.S  <- rownames(S)
        ID.Vr <- rownames(V)
        ID.Vc <- colnames(V)
        
        if (length(setdiff(ID, ID.S)) != 0)
            return(print("ERROR: some ID genotype are not found in ID structure"))
        
        if (length(setdiff(ID, ID.Vr)) != 0 | length(setdiff(ID, ID.Vc)) != 0 )
            return(print("ERROR: some ID genotype are not found in ID genotypic variance-covariance"))
        
        alignS            <- merge(donnees, S, by.x = 0, by.y = 0, all.x = TRUE, sort = FALSE)
        donnees           <- alignS[, 2:(ncol(donnees) + 1)]
        rownames(donnees) <- alignS[, 1]
        S                 <- as.data.frame(alignS[, (ncol(donnees) + 2):ncol(alignS)])
        rownames(S)       <- alignS[, 1]
        colnames(S)       <- colnames(alignS)[(ncol(donnees)+2):ncol(alignS)]
        alignV            <- merge(V, donnees, by.x = 0, by.y = 0, sort = FALSE)
        donnees           <- alignV[,(ncol(V)+2):ncol(alignV)]
        rownames(donnees) <- alignV[,1]
        rownames(alignV)  <- alignV[,1]
        V                 <- alignV[,-1]
        V                 <- t(V)
        V                 <- merge(V, donnees, by.x = 0, by.y = 0, sort = FALSE)
        rownames(V)       <- V[, 1]
        V                 <- V[, 2:(nrow(V) + 1)]
        V                 <- as.matrix(V)
        S                 <- as.matrix(S)
        
    } else {
        if (!S.yes) {
            ID   <- rownames(donnees)
            ID.S <- rownames(S)
            
            if (length(setdiff(ID,ID.S)) != 0)
                return(print("ERROR: some ID genotype are not found in ID structure"))
            
            alignS            <- merge(donnees, S, by.x = 0, by.y = 0, all.x = TRUE, sort = FALSE)
            donnees           <- alignS[, 2:(ncol(donnees) + 1)]
            rownames(donnees) <- alignS[, 1]
            S                 <- as.data.frame(alignS[, (ncol(donnees) + 2):ncol(alignS)])
            rownames(S)       <- alignS[, 1]
            colnames(S)       <- colnames(alignS)[(ncol(donnees) + 2):ncol(alignS)]
            S                 <- as.matrix(S)
        }
        
        if (!V.yes) {
            ID    <- rownames(donnees)
            ID.Vr <- rownames(V)
            ID.Vc <- colnames(V)
            
            if (length(setdiff(ID, ID.Vr)) != 0 | length(setdiff(ID, ID.Vc)) != 0 )
                return(print("ERROR: some ID genotype are not found in ID genotypic variance-covariance"))
            
            alignV            <- merge(V, donnees, by.x = 0, by.y = 0, sort = FALSE)
            donnees           <- alignV[, (ncol(V) + 2):ncol(alignV)]
            rownames(donnees) <- alignV[, 1]
            rownames(alignV)  <- alignV[, 1]
            V                 <- alignV[, -1]
            V                 <- t(V)
            V                 <- merge(V, donnees, by.x = 0, by.y = 0, sort = FALSE)
            rownames(V)       <- V[, 1]
            V                 <- V[, 2:(nrow(V) + 1)]
            V                 <- as.matrix(V)
        }
    }
    
    donnees <- as.matrix(donnees)
    
    if (supinfo)
        info <- apply(donnees, data = data, 2, Info.Locus)
    
    
    if (!na.presence & !V.yes) {
        V_inv           <- Inv.proj.matrix.sdp(V) 
        rownames(V_inv) <- rownames(V)
        colnames(V_inv) <- colnames(V)
    }
    
    RES <- mclapply(1:(ncol(donnees) - 1), function(i) {
        res <- sapply((i + 1):(ncol(donnees)), function(j) {
            li   <- donnees[, i]
            lj   <- donnees[, j]
            lij  <- cbind(li, lj)
            loc1 <- c(loc1, mark[i])
            loc2 <- c(loc2, mark[j])
            
            # options(warn = -1)
            M.r2 <- Measure.R2(biloci = lij, na.presence = na.presence)
            # options(warn = 0)
            
            if (!S.yes)
                M.r2s <- Measure.R2S(biloci = lij, struc = S, na.presence = na.presence)
            
            if (!V.yes) {
                if(!na.presence)
                    M.r2v <- Measure.R2V(biloci = lij, V = V, na.presence = na.presence, V_inv = V_inv)
                else
                    M.r2v <- Measure.R2V(biloci = lij, V = V, na.presence = na.presence)
            }
            
            if (!S.yes & !V.yes) {
                if(!na.presence)
                    M.r2vs <- Measure.R2VS(biloci = lij, V = V, struc = S, na.presence = na.presence, V_inv = V_inv)
                else
                    M.r2vs <- Measure.R2VS(biloci = lij, V = V, struc = S, na.presence = na.presence)
            }
            
            r2   <- c(r2, M.r2)
            r2s  <- c(r2s, M.r2s)
            r2v  <- c(r2v, M.r2v)
            r2vs <- c(r2vs, M.r2vs)
            
            if (supinfo) {
                M.MAF.loc1        <- NA
                M.heterofreq.loc1 <- NA
                M.NAfreq.loc1     <- NA
                M.MAF.loc2        <- NA
                M.heterofreq.loc2 <- NA
                M.NAfreq.loc2     <- NA
                
                supinfo.li        <- info[, which(colnames(info) == mark[i])]
                supinfo.lj        <- info[, which(colnames(info) == mark[j])]
                M.MAF.loc1        <- supinfo.li[1]
                M.heterofreq.loc1 <- supinfo.li[2]
                M.NAfreq.loc1     <- supinfo.li[3]
                M.MAF.loc2        <- supinfo.lj[1]
                M.heterofreq.loc2 <- supinfo.lj[2]
                M.NAfreq.loc2     <- supinfo.lj[3]
                
                MAF.loc1        <- c(MAF.loc1, M.MAF.loc1)
                heterofreq.loc1 <- c(heterofreq.loc1, M.heterofreq.loc1)
                NAfreq.loc1     <- c(NAfreq.loc1, M.NAfreq.loc1)
                MAF.loc2        <- c(MAF.loc2, M.MAF.loc2)
                heterofreq.loc2 <- c(heterofreq.loc2, M.heterofreq.loc2)
                NAfreq.loc2     <- c(NAfreq.loc2, M.NAfreq.loc2)
            }
            
            res <- c(as.matrix(loc1), as.matrix(loc2), r2, r2s, r2v, r2vs,
                     MAF.loc1, heterofreq.loc1, NAfreq.loc1,
                     MAF.loc2, heterofreq.loc2, NAfreq.loc2)
        })
        
        RES <- cbind(RES, res)
    })
    
    RES <- matrix(unlist(RES), nrow = length(RES[[1]][, 1]))
    
    if (!S.yes & !V.yes) {
        result <- data.frame(RES[1, ], RES[2, ], as.numeric(RES[3, ]), as.numeric(RES[5, ]), as.numeric(RES[4, ]), as.numeric(RES[6, ]))
        colnames(result) <- c("loc1", "loc2", "r2", "r2v", "r2s", "r2vs")
        
    } else {
        if (!S.yes) {
            result <- data.frame(RES[1, ], RES[2, ], as.numeric(RES[3, ]), as.numeric(RES[4, ]))
            colnames(result) <- c("loc1","loc2", "r2", "r2s")
        } else {
            if (!V.yes) {
                result <- data.frame(RES[1, ], RES[2, ], as.numeric(RES[3, ]), as.numeric(RES[5, ]))
                colnames(result) <- c("loc1", "loc2", "r2", "r2v")
            } else {
                result <- data.frame(RES[1, ], RES[2, ], as.numeric(RES[3, ]))
                colnames(result) <- c("loc1", "loc2", "r2")
            }
        }
    }
    
    if (supinfo) {
        Info <- data.frame(as.numeric(RES[7, ]), as.numeric(RES[8, ]), as.numeric(RES[9, ]), as.numeric(RES[10, ]), as.numeric(RES[11, ]), as.numeric(RES[12, ]))
        result <- as.data.frame(c(result, Info))
        l <- length(colnames(result))
        colnames(result)[(l - 5):(l)] <- c("MAF.loc1", "heterofreq.loc1", "NAfreq.loc1",
                                           "MAF.loc2", "heterofreq.loc2", "NAfreq.loc2")
        
    }	
    result <- as.data.frame(result)
    result
}
