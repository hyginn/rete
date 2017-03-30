# diffuseTools.R
#
# Utility functions for heat diffusion

.nonEqSim <- function(vTable, eTable, N, q, seed, silent) {
    # Purpose:
    #     ...
    #
    # Parameters:
    #     vTable: vertex information
    #     eTable: edge information
    #     N:      ...
    #     q:      ...
    #     seed:   ...
    #     silent: ...
    # Details:
    #     ....
    #

    # Value:
    #     gList: updated eTable and vTable

    # ==== SETUP METADATA ======================================================

    set.seed(seed)
    nE <- nrow(eTable)
    N <- N * nE

    # For randomly chosen edges of G repeat N times:


    # Compute the difference in heat between the incident vertices.
    # Take a fraction q of heat from one vertex and add it to the other.
    # Record the amount of heat that has been passed as an attribute of the
    #  edge in the data structure that stores "flux".

    for (i in 1:N) { # run for N steps ...
        if (!silent) { .pBar(i, N) }

        thisE <- sample(1:nE, 1)   # pick a random edge
        vT <- eTable$a[thisE]      # tail vertex
        vH <- eTable$b[thisE]      # head Vertex
        hvT <- vTable$heat[vT]     # heat of Tail Vertex
        hvH <- vTable$heat[vH]     # heat of Head Vertex
        if (hvT > hvH) {           # ... equilibrate
            dH <- q * eTable$W[thisE] * ((hvT - hvH)/2)   # heat difference
                                                          # portion:
            vTable$heat[vT] <- vTable$heat[vT] - dH  # flow, from here ...
            vTable$heat[vH] <- vTable$heat[vH] + dH  # ... to there,
            eTable$flux[thisE] <- eTable$flux[thisE] + dH # add it to the flux.
        }
    }

    return(list(eTable = eTable, vTable = vTable))
}

# [END]
