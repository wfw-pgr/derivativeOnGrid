import nkUtilities.generate__sampleProfile as gsp

# ========================================================= #
# ===  generate Sample                                  === #
# ========================================================= #
def genSample():
    outFile     = "dat/field.dat"
    x1MinMaxNum = [ 0.0, 1.0, 21 ]
    x2MinMaxNum = [ 0.0, 1.0, 21 ]
    x3MinMaxNum = [ 0.0, 1.0, 21 ]
    ret = gsp.generate__sampleProfile( x1MinMaxNum=x1MinMaxNum, x2MinMaxNum=x2MinMaxNum, \
                                       x3MinMaxNum=x3MinMaxNum, outFile    =outFile      )


# ======================================== #
# ===  実行部                          === #
# ======================================== #
if ( __name__=="__main__" ):
    genSample()
