getXKX_hdf5 <- function(filename, XX, K, inv, lambda, mc.cores=1){
  # ans <- apply(1:length(XX), getXKX_bd.i, filename=filename, XX=XX, K=K, inv=inv,
  #                 lambda=lambda, nthreads=mc.cores)
  # ans

  ######
  ######  FALTA REPASSAR TOTA LA FUNCIÓ I ACABAR-HO DE PROGRAMAR PER CALCULAR
  ######  EL PRODUCTE SCALAR + LA INVERSIÓ CHOLEVSKY !!! (ELSE)
  ######

  # bdapply_Function_hdf5(filename = filename,
  #                       group = "X",datasets = XX,
  #                       b_group = "K",b_datasets = K,
  #                       outgroup = "XK",func = "blockmult",
  #                       transp_bdataset = TRUE,
  #                       force = T)
  #
  # XK <- bdgetDatasetsList_hdf5(filename = filename, group = "XK")
  # bdapply_Function_hdf5(filename = filename,
  #                       group = "XK",datasets = XK,
  #                       b_group = "X",b_datasets = XX,
  #                       outgroup = "M",func = "blockmult",
  #                       transp_bdataset = TRUE,
  #                       force = T)


    bdapply_Function_hdf5(filename = filename,
                          group = "X",datasets = XX,
                          outgroup = "M",func = "tCrossProd",
                          force = T)

    M <- bdgetDatasetsList_hdf5(filename = filename, group = "M")

    if (inv==1) { # solve

    bdapply_Function_hdf5(filename = filename,
                          group = "M",datasets = M,
                          outgroup = "XKX",func = "invChol",
                          force = T)

    } else if (inv==2) { # penalized


    #. A programar .# xkx <- bdInvCholesky(M + bdScalarwproduct(diag(nrow(M)), lambda[i], "wX"))

      # ==>
      tmp <- bdgetDatasetsList_hdf5(filename = filename, group = "M")

      #
      # Fer el pas inicial i emmagatzemar al grup tmp
      #
      sapply(1:length(tmp), function(i) {
                  #..# tmpResult <- bdScalarwproduct( bdgetDiagonal_hdf5(filename, "M", tmp[i]) , lambda[i], "wX")
                dims <- BigDataStatMeth::bdgetDim_hdf5(filename, paste0("M/",tmp[i]))
                tmpResult <- bdScalarwproduct( diag(dims[1]) , lambda[i], "wX")
          browser()
                  # aa <- tmp[i] + bdScalarwproduct( bdgetDiagonal_hdf5(filename, "M", tmp[i]) , lambda[i], "wX")
                  bdAdd_hdf5_matrix( tmpResult, filename, "tmp", paste0(tmp[i],"scalarx"), force = T)

                  print(paste0("Sumatori de ", tmp[i]))

                  bdblockSum_hdf5( filename = filename,
                                   group = "M", a = tmp[i],
                                   groupB = "tmp", b = paste0(tmp[i],"scalarx"),
                                   outgroup = "tmp", outdataset = tmp[i],
                                   block_size = 7500 )

          } )

      print("Just ABANS del aplay de cholesky")
      print("Variables valen : ")
      print("outgroup = tmp")
      print(paste0("datasets = ", tmp, " -- "))
      print(paste0("filename = ", filename))
      browser()
      bdapply_Function_hdf5(filename = filename,
                            group = "tmp",datasets = tmp,
                            outgroup = "XKX",func = "invChol",
                            force = T)

      print("Just DESPRES del aplay de cholesky")
      # <==

    } else {
    stop("need correct method")
    }

    # if (inv==1) # solve
    #   xkx <- bdInvCholesky(M)
    # else if (inv==2) # penalized
    #   xkx <- bdInvCholesky(M + bdScalarwproduct(diag(nrow(M)), lambda[i], "wX"))
    # else
    #   stop("need correct method")

}

getXKX_hdf5.i <- function(filename, i, XX, K, inv, lambda, nthreads)
{

  #.Ok.# xk <- BigDataStatMeth::bdblockmult(t(XX[[i]]),K[[i]], onmemory = T, threads = nthreads)
  #.Ok.# M <- BigDataStatMeth::bdblockmult(xk,XX[[i]], onmemory = T, threads = nthreads)

  if (inv==1) # solve
    xkx <- bdInvCholesky(M)
  else if (inv==2) # penalized
    xkx <- bdInvCholesky(M + bdScalarwproduct(diag(nrow(M)), lambda[i], "wX"))
  else
    stop("need correct method")

  ans <- list(xkx=xkx, xk=xk)
  ans
}




## El que faltaria fer es, al realitzar el sumatori,  mirar si el que em de llegir
## es una matriu o un vector.
##  - Si es un vector, replicar el vector per aconseguir la mida de la matriu
##  - Si es una matriu, doncs res, llegir la matriu i llestos com fins ara.

