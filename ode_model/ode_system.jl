function ODEfunc!(dydt, y, p, t)
    #Initiation
    scaling_fact, scaling_fact2 = p
    dydt[1] =(kf_C3h2o * y[88] - kf_C3h2oB * y[1] * y[90] + kb_C3h2oB * y[2]
    +kb_C3bH * y[46] - kf_C3bH * y[1] * y[95] - kf_C3bCR1 * y[1] * y[96]
    +kb_C3bCR1 * y[47] + kb_C3h2oBb * y[3] + kb_C3bBbH * y[48])

    dydt[2] = (kf_C3h2oB * y[1] * y[90] - kb_C3h2oB * y[2] - 
                kcat_C3h2oB * y[91] * y[2]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21]) )

    dydt[3] = (-kf_C3bH * y[3] * y[95] - kb_C3h2oBb * y[3] + 
                kcat_C3h2oB * y[91] * y[2]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21]) )
                 
    dydt[4] = (kf_fC3b * y[7] * y[8] + kf_fC3b * y[13] * y[8] - kf_C3Bb * y[4] * y[90] + kb_C3Bb * y[5] - 
    kf_C3bH * y[4] * y[95] + kb_C3bH * y[62] - kf_C3bP * y[4] * y[12] + kb_C3bP * y[19] - 
    kf_C3bCR1 * y[4] * y[96] + kb_C3bCR1 * y[64] + kb_C3bBb * y[6] + kb_C3bBbH * y[67] + kb_C3bBbCR1 * y[68] - 
    kf_C3bBbC3b * y[53] * y[4] - kf_C3bBbC3b * y[4] * y[15] - kf_C3bBbC3b * y[4] * y[18] - kf_C3bBbC3b * y[4] * y[21] + 
    kb_C3bBbC3b * y[69] + kb_C3bBbC3b * y[27] + kb_C3bBbCR1 * y[71] + kb_C3bBbDAF * y[72] + 
    kb_C3bBbH * y[70] + kb_p⁺su * y[19] + kb_p⁺su * y[20] + kb_p⁺su * y[21] +
    kb_p⁺su * y[32] + kb_p⁺su * y[36] + kb_p⁺su * y[38] + kb_p⁺su * y[34])
     
    dydt[5] = (kf_C3Bb * y[4] * y[90]  - kb_C3Bb * y[5] -
    kcat_C3h2oB * y[91] * y[5]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21]) )      

    dydt[6] = ( -kb_C3bBb * y[6] - kf_C3bH * y[6] * y[95] -  kf_C3bCR1 * y[6] * y[96] + 
    kcat_C3h2oB * y[91] * y[5]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21]) )   

    dydt[7] = (kcat_C3h2oBb * y[88] * y[3]/(km_C3h2oBb + y[88])  + kcat_C3bBb * y[88] * y[6]/(km_C3bBb + y[88]) - kf_fC3b * y[7] * y[8]
     -kf_C3bP * y[7] * y[12] - kf_C3Bsu * y[7] * y[50]  - kf_C3Bsu * y[7] * y[10])
    
    dydt[8] = -kf_fC3b * y[7] * y[8] -kf_fC3b * y[13] * y[8] -kf_fC3b * y[49] * y[8] 

    #Amplification on pathogen
    dydt[9] =  kf_p⁺re * y[9]

    dydt[10] = (-kf_C3Bsu * y[7] * y[10] - kf_pC3b * y[10] * y[13] * scaling_fact - 
                    kf_p⁺su * y[94] * y[10] - kf_C5b7su * y[39] * y[10] * scaling_fact2)


    dydt[11]= ( kf_C3Bsu * y[7] * y[10] + kf_pC3b * y[10] * y[13] * scaling_fact - 
            kf_C3Bb * y[11] * y[90] + kb_C3Bb * y[14] - kf_C3bP * y[11] * y[94] + kb_C3bP * y[19] +
            kb_C3bBb * y[15]   + kb_C3bBbP * y[18] - kf_C3bCR1 * y[11] * y[96] +
            kb_C3bCR1 * y[22])
    
    dydt[12] = (kf_p⁺su * y[94] * y[10]   - kb_p⁺su * y[12] - y[12] * kf_C3bP * (y[13] + y[7] + y[4]) +
    kb_C3bP * (y[19] + y[32] + y[20] + y[21] + y[34]+ y[36]+ y[38]))

    dydt[13] =  (kcat_C3bBb * y[88] * y[15]/(km_C3bBb+ y[88]) + kcat_C3bBbP * y[88] * y[18]/(km_C3bBbP+ y[88]) + 
    kcat_C3bBbP * y[88] * y[21]/(km_C3bBbP+ y[88]) -
    kf_fC3b * y[13] * y[8] - kf_pC3b * y[10] * y[13] * scaling_fact -  y[12] * kf_C3bP * y[13] + kb_C3bP * y[19] - 
    kf_C3bBbC3b * y[13] * y[15] - kf_C3bBbC3b * y[13] * y[18] - kf_C3bBbC3b * y[13] * y[21]) 

    dydt[14] = (kf_C3Bb * y[11] * y[90] - kb_C3Bb * y[14]  - kf_C3bP * y[14] * y[94] + kb_C3bP * y[20]  - 
    kcat_C3h2oB * y[91] * y[14]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21]) )  


    dydt[15] = (-kf_C3bP * y[15] * y[94] -kf_C3bP * y[15] * y[93]  + kb_C3bBbP * y[18] - kf_C3bBbC3b * y[4] * y[15] -
    kf_C3bBbC3b * y[13] * y[15] + kb_C3bBbC3b * y[27] +
    kcat_C3h2oB * y[91] * y[14]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21]) )  

    dydt[16] = (kf_C3bP * y[11] * y[94] + kb_C3bP * y[20])

    dydt[17] = (kf_C3bP * y[19] * y[90] + kf_C3bP * y[14] * y[94] - kb_C3bP * y[20] - 
    kcat_C3h2oB * y[91] * y[17]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21]))

    dydt[18] = (kf_C3bP * y[15] * y[94] + kf_C3bP * y[15] * y[93] - kb_C3bBbP * y[18] - 
    kf_C3bBbC3b * y[4] * y[18] - kf_C3bBbC3b * y[13] * y[18] + 
    kcat_C3h2oB * y[91] * y[17]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21]))

    dydt[19] = (kf_C3bP * y[12] * (y[13] + y[7] + y[4]) - kb_C3bP * y[19] - kb_p⁺su * y[19] - kf_C3bB * y[19] * y[90])

    dydt[20] = (kf_C3bB * y[19] * y[90] - kb_C3bP * y[20] - kb_p⁺su * y[20] - 
    kcat_C3h2oB * y[91] * y[20]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21]))

    dydt[21] = (-kb_C3bBbP * y[21] - kf_C3bBbC3b * y[4] * y[21] - kf_C3bBbC3b * y[13] * y[21] -
    kb_p⁺su * y[21] + 
    kcat_C3h2oB * y[91] * y[20]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21]))

    dydt[22] = (kf_C3bCR1 * y[11] * y[96] - kb_C3bCR1 * y[22] - 
    kcat_C3bH * y[92] * y[22]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]))

    dydt[23] = (-kf_iC3bCR1 * y[23] * y[96] + kb_iC3bCR1 * y[24] - kf_iC3bP * y[23] * y[93] + kb_iC3bP * y[26] + 
    kcat_C3bH * y[92] * y[22]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]))

    dydt[24] = (kf_iC3bCR1 * y[23] * y[96] - kb_iC3bCR1 * y[24] -
    kcat_C3bH * y[92] * y[24]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]))

    dydt[25] = ( kcat_C3bH * y[92] * y[24]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]))

    dydt[26] = (-kb_iC3bP * y[26] + kf_iC3bP * y[23] * y[93])

    #Termination Reactions
    dydt[27] = (kf_C3bBbC3b * y[4] * y[18] + kf_C3bBbC3b * y[13] * y[18] - kb_C3bBbC3b * y[27] - 
    kf_C3bP * y[27] * (y[93] + y[94]) - kf_C3bBbC3bC5 * y[27] * y[98] + kb_C3bBbC3bC5 * y[28] + 
    kb_C5b * y[29] + kb_C3bP * y[31] + kf_C5b7 * y[30] * y[102]) 

    dydt[28] = (kf_C3bBbC3bC5 * y[27] * y[98] - kb_C3bBbC3bC5 * y[28] + kb_C3bP * y[33] - 
    kcat_C3bBbC3b * y[28]) 

    dydt[29] = (kcat_C3bBbC3b * y[28] - kb_C5b * y[29] + kb_C3bP * y[35] - 
    kf_C3bBbC3bC5bC6 * y[29] * y[101] + kb_C3bBbC3bC5bC6 * y[30])                                             


    dydt[30] = (kf_C3bBbC3bC5bC6 * y[29] * y[101] - kb_C3bBbC3bC5bC6 * y[30] +
    kb_C3bP * y[37] - kf_C5b7 * y[30] * y[102]) 

    dydt[31] = (kf_C3bBbC3b * (y[4] + y[13]) * y[18] + kf_C3bP * y[27] * (y[93] + y[94])
    + kf_C5b7 * y[30] * y[102] - kb_C3bP * y[31] + kb_C5b * y[35] - 
    kf_C3bBbC3bC5 * y[31] * y[98] + kb_C3bBbC3bC5 * y[33]) 

    dydt[32] = (kf_C3bBbC3b * (y[4] + y[13]) * y[21] - kb_C3bP * y[32] + 
    kf_C5b7 * y[38] * y[102] - kf_C3bBbC3bC5 * y[32] * y[98] + kb_C3bBbC3bC5 * y[34] + 
    kb_C5b * y[36] - kb_p⁺su * y[32]) 


    dydt[33] = (kf_C3bBbC3bC5 * y[31] * y[98] - kb_C3bBbC3bC5 * y[33] - 
    kcat_C3bBbC3b * y[33] - kb_C3bP * y[33]) 

    dydt[34] = (kf_C3bBbC3bC5 * y[32] * y[98] - kb_C3bBbC3bC5 * y[34] - 
    kcat_C3bBbC3b * y[34] - kb_p⁺su * y[34]) 

    dydt[35] = (kf_C3bBbC3bC5 * y[31] * y[98] - kb_C5b * y[35] - kb_C3bP * y[35] - 
    kf_C3bBbC3bC5bC6 * y[35] * y[101] + kb_C3bBbC3bC5bC6 * y[37]) 

    dydt[36] = (kf_C3bBbC3bC5 * y[32] * y[98] - kb_C5b * y[36] - 
    kf_C3bBbC3bC5bC6 * y[36] * y[101] + kb_C3bBbC3bC5bC6 * y[38] - 
    kb_p⁺su * y[36]) 

    dydt[37] = (kf_C3bBbC3bC5bC6 * y[35] * y[101] - kb_C3bBbC3bC5bC6 * y[37] - 
    kf_C5b7 * y[37] * y[102] - kb_C3bP * y[37]) 

    dydt[38] = (kf_C3bBbC3bC5bC6 * y[36] * y[101] - kb_C3bBbC3bC5bC6 * y[38] - 
    kf_C5b7 * y[38] * y[102] - kb_p⁺su * y[38]) 

    dydt[39] = (kf_C5b7 * y[102] * (y[38] + y[37] + y[30]) - 
    kf_C5b7su * y[39] * y[10] * scaling_fact2 - kf_mi * y[39] - 
    kf_C5b8 * y[39] * y[103] + kb_C5b8 * y[41] - kb_C5b7 * y[39] -
    kf_CnC5b7 * y[39] * y[106]  + kb_CnC5b7 * y[85] - kf_VnC5b7 * y[39] * y[105] + kb_VnC5b7 * y[82]) 

    dydt[40] = (kf_mi * (y[39] + y[76])) 

    dydt[41] = (kf_C5b8 * y[103] * (y[76] + y[39]) - kb_C5b8 * y[41] - kf_C5b9 * y[41] * y[104] + kb_C5b9 * y[42]) 

    dydt[42] = ( kf_C5b9 * y[41] * y[104] - kb_C5b9 * y[42]) 

    dydt[43] = (kf_C5b7su * y[39] * y[10] * scaling_fact2 - kf_C5b8 * y[43] * y[103]) 

    dydt[44] = (kf_C5b8 * y[43] * y[103] - kf_C5b9 * y[44] * y[104]) 

    dydt[45] = (kf_C5b9 * y[44] * y[104])  
    
    # Regulation terms
    dydt[46] = kf_C3bH * y[1] * y[95] - kb_C3bH * y[46]

    dydt[47] = kf_C3bCR1 * y[1] * y[96] - kb_C3bCR1 * y[47]
    
    dydt[48] = kf_C3bH * y[3] * y[95] - kb_C3bBbH * y[48]
    
    dydt[49] = (kcat_C3bBb * y[88] * y[53]/(km_C3bBb + y[88]) - kf_fC3b * y[49] * y[8] - 
    kf_C3Bsu * y[49] * y[50] * 15.6 - kf_C3bBbC3b * y[53] * y[49])
    
    dydt[50] = (-kf_C3Bsu * y[7] * y[50]  - kf_C3Bsu * y[49] * y[50] * 15.6 - kf_C5b7su * y[76] * y[50] * 4*1.4)
    
    dydt[51] = (kf_C3Bsu * y[7] * y[50] + kf_C3Bsu * y[49]  * 15.6 * y[50] - kf_C3bCR1 * y[51] * y[96] + 
    kb_C3bCR1 * y[57] - kf_C3bH * y[51] * y[95] + kb_C3bH_ho * y[59] - kf_C3bB * y[51] * y[90] + kb_C3bB * y[52] + 
    kb_C3bBbH * y[54] + kb_C3bBbCR1 * y[55] + kb_C3bBbDAF * y[56] + kb_C3bBbH * y[70] + 
    kb_C3bBbC3b * y[69] + kb_C3bBbCR1 * y[71] + kb_C3bBbDAF * y[72])
    
    dydt[52] = (kf_C3bB * y[51] * y[90] - kb_C3bB * y[52] - 
    kcat_C3h2oB * y[91] * y[52]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21]) )
    
    dydt[53] = (- kf_C3bBbC3b * y[53] * (y[4] + y[49]) - kf_C3bH * y[53] * y[95] - kf_C3bCR1 * y[53] * y[96] - 
    kf_C3bBbDAF * y[53] * y[97] + 
    kcat_C3h2oB * y[91] * y[52]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21]) )
    
    dydt[54] = (kf_C3bH * y[53] * y[95] - kb_C3bBbH * y[54])

    dydt[55] = (kf_C3bCR1 * y[53] * y[96] - kb_C3bBbCR1 * y[55])
    
    dydt[56] = (kf_C3bBbDAF * y[53] * y[97] - kb_C3bBbDAF * y[56])
    
    dydt[57] = (kf_C3bCR1 * y[51] * y[96] - kb_C3bCR1 * y[57] - 
    kcat_C3bH * y[92] * y[57]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]))
    
    dydt[58] = (-kf_iC3bCR1 * y[58] * y[96] + kb_iC3bCR1 * y[60] + 
    kcat_C3bH * y[92] * y[59]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]) + 
    kcat_C3bH * y[92] * y[57]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]))

    dydt[59] = (kf_C3bH * y[51] * y[95] - kb_C3bH_ho * y[59] - 
    kcat_C3bH * y[92] * y[59]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]) )

    dydt[60] = (kf_iC3bCR1 * y[58] * y[96] - kb_iC3bCR1 * y[60] - 
    kcat_C3bH * y[92] * y[60]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]))

    dydt[61] = ( kcat_C3bH * y[92] * y[60]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]))
 
    dydt[62] = ( kf_C3bH * y[4] * y[95] - kb_C3bH * y[62] - 
    kcat_C3bH * y[92] * y[62]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]))

    dydt[63] = (-kf_iC3bCR1 * y[63] * y[96] + kb_iC3bCR1 * y[65] + 
    kcat_C3bH * y[92] * y[62]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]) + 
    kcat_C3bH * y[92] * y[64]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]))
    
    dydt[64] = (kf_C3bCR1 * y[4] * y[96] - kb_C3bCR1 * y[64] - 
    kcat_C3bH * y[92] * y[64]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]))
    
    dydt[65] = (kf_iC3bCR1 * y[63] * y[96] - kb_iC3bCR1 * y[65] - 
    kcat_C3bH * y[92] * y[65]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]))
    
    dydt[66] =  ( kcat_C3bH * y[92] * y[65]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]))
    
    dydt[67] = (kf_C3bH * y[6] * y[95] - kb_C3bBbH * y[67])
    
    dydt[68] = (kf_C3bCR1 * y[6] * y[96] - kb_C3bBbCR1 * y[68])
    
    dydt[69] = (kf_C3bBbC3b * y[49] * y[53] + kf_C3bBbC3b * y[4] * y[53] - kb_C3bBbC3b * y[69] - 
    kf_C3bCR1 * y[69] * y[96]  - kf_C3bBbDAF * y[69] * y[97] - kf_C3bH * y[69] * y[95] - 
    kf_C3bBbC3bC5 * y[69] * y[98] + kb_C3bBbC3bC5 * y[73] + kb_C5b * y[74] + kf_C5b7 * y[75] * y[102]) 

    dydt[70] = (kf_C3bH * y[69] * y[95] - kb_C3bBbH * y[70])
    
    dydt[71] = (kf_C3bCR1 * y[69] * y[96] - kb_C3bBbCR1 * y[71])
    
    dydt[72] = (kf_C3bBbDAF * y[69] * y[97] - kb_C3bBbDAF * y[72])
    
    dydt[73] = (kf_C3bBbC3bC5 * y[69] * y[98] - kb_C3bBbC3bC5 * y[73]  - 
    kcat_C3bBbC3b * y[73]) 
    
    dydt[74] = (kcat_C3bBbC3b * y[73] - kb_C5b * y[74] - 
    kf_C3bBbC3bC5bC6 * y[74] * y[101] + kb_C3bBbC3bC5bC6 * y[75])
    
    dydt[75] = (kf_C3bBbC3bC5bC6 * y[74] * y[101] - kb_C3bBbC3bC5bC6 * y[75] - 
    kf_C5b7 * y[75] * y[102]) 
    
    dydt[76] = (kf_C5b7 * y[75] * y[102] - kb_C5b7 * y[76] - kf_C5b7su * y[76]  * y[50] * 4*1.4 - (kf_mi * y[76]) -
    kf_C5b8 * y[76] * y[103] - kf_CnC5b7 * y[76] * y[106] + kb_CnC5b7 * y[85] - kf_VnC5b7 * y[76] * y[105] + 
    kb_VnC5b7 * y[82] )
    
    dydt[77] = (kf_C5b7su * y[76]  * y[50] * 4*1.4 - kf_C5b8 * y[77] * y[103])
    
    dydt[78] = (kf_C5b8 * y[77] * y[103] - kf_C5b9 * y[78] * y[104])
    
    dydt[79] = (kf_C5b9 * y[78] * y[104] - kf_C5b9 * y[79] * y[104] - kf_CD59C5b9 * y[79] * y[107] + kb_CD59C5b9 * y[81])
    
    dydt[80] = (kf_C5b9 * y[79] * y[104])
    
    dydt[81] = (kf_CD59C5b9 * y[79] * y[107] - kb_CD59C5b9 * y[81])

    dydt[82] = (kf_VnC5b7 * y[76] * y[105] - kb_VnC5b7 * y[82] + kf_VnC5b7 * y[39] * y[105] - kf_VnC5b8 * y[82] * y[103] +
    kb_VnC5b8 * y[83])
    
    dydt[83] = (kf_VnC5b8 * y[82] * y[103] - kb_VnC5b8 * y[83] - kf_VnC5b9 * y[83] * y[104] + kb_VnC5b9 * y[84])
    
    dydt[84] = (kf_VnC5b9 * y[83] * y[104] - kb_VnC5b9 * y[84])
    
    dydt[85] = (kf_CnC5b7 * y[76] * y[106] - kb_CnC5b7 * y[85] + kf_CnC5b7 * y[39] * y[106] - kf_CnC5b8 * y[85] * y[103] +
    kb_CnC5b8 * y[86])
    
    dydt[86] = (kf_CnC5b8 * y[85] * y[103] - kb_CnC5b8 * y[86] - kf_CnC5b9 * y[86] * y[104] + kb_CnC5b9 * y[87])
    
    dydt[87] = (kf_CnC5b9 * y[86] * y[104] - kb_CnC5b9 * y[87])

    # Complement Proteins (Host cell and fluid state)
    dydt[88] = (-kf_C3h2o * y[88] - kcat_C3h2oBb * y[88] * y[3]/(km_C3h2oBb + y[88])  - kcat_C3bBb * y[88] * y[6]/(km_C3bBb + y[88]) -
    kcat_C3bBb * y[88] * y[15]/(km_C3bBb+ y[88]) - kcat_C3bBbP * y[88] * y[18]/(km_C3bBbP+ y[88]) - 
    kcat_C3bBbP * y[88] * y[21]/(km_C3bBbP+ y[88]) - kcat_C3bBb * y[88] * y[53]/(km_C3bBb + y[88]) )  
    
    dydt[89] = (kcat_C3h2oBb * y[88] * y[3]/(km_C3h2oBb + y[88]) + kcat_C3bBb * y[88] * y[6]/(km_C3bBb + y[88]) +
    kcat_C3bBb * y[88] * y[15]/(km_C3bBb+ y[88]) + kcat_C3bBbP * y[88] * y[18]/(km_C3bBbP+ y[88]) + 
    kcat_C3bBbP * y[88] * y[21]/(km_C3bBbP+ y[88]) + kcat_C3bBb * y[88] * y[53]/(km_C3bBb + y[88])) 

    # dydt[89] = (kcat_C3bBbP * y[88] * y[18]/(km_C3bBbP+ y[88])) 

    dydt[90] = (-kf_C3h2oB * y[1] * y[90] + kb_C3h2oB * y[2] - kf_C3Bb * y[4] * y[90] + kb_C3Bb * y[5] - 
    kf_C3Bb * y[11] * y[90] + kb_C3Bb * y[14] - kf_C3bB * y[19] * y[90] - kf_C3bB * y[51] * y[90] + kb_C3bB * y[52])  

    dydt[91] = (- kcat_C3h2oB * y[91] * y[2]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21]) -
    kcat_C3h2oB * y[91] * y[5]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21]) -
    kcat_C3h2oB * y[91] * y[52]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21]) - 
    kcat_C3h2oB * y[91] * y[14]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21]) - 
    kcat_C3h2oB * y[91] * y[20]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21]) - 
    kcat_C3h2oB * y[91] * y[17]/(km_C3h2oB + y[5] + y[2] + y[52] + y[14] + y[17] + y[21])) 

    dydt[92] = (-kcat_C3bH * y[92] * y[22]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24])-
    kcat_C3bH * y[92] * y[24]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]) - 
    kcat_C3bH * y[92] * y[57]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]) - 
    kcat_C3bH * y[92] * y[59]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]) - 
    kcat_C3bH * y[92] * y[60]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]) - 
    kcat_C3bH * y[92] * y[62]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]) - 
    kcat_C3bH * y[92] * y[64]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]) - 
    kcat_C3bH * y[92] * y[65]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24])) 

    dydt[93] = (-kf_C3bP * y[15] * y[93] + kb_C3bBbP * y[18] - kf_C3bP * y[27] * y[93] - kf_iC3bP * y[23] * y[93] +
    kb_iC3bP * y[26] + kb_C3bP * y[31] )  

    dydt[94] = (kf_p⁺re * y[9] - kf_p⁺su * y[94] * y[10] + kb_p⁺su * y[12] - kf_C3bP * y[11] * y[94] + kb_C3bP * y[19] - 
    kf_C3bP * y[14] * y[94] + kb_C3bP * y[20] -kf_C3bP * y[15] * y[94] + kb_C3bBbP * y[21] - 
    kf_C3bP * y[27] * y[94] + kb_C3bP * y[32] + kb_p⁺su * y[19] + kb_p⁺su * y[20] + kb_p⁺su * y[21] +
    kb_p⁺su * y[32] + kb_p⁺su * y[36] + kb_p⁺su * y[38] + kb_p⁺su * y[34]) 


    dydt[95] = (-kf_C3bH * y[1] * y[95] + kb_C3bH * y[46] - kf_C3bH * y[3] * y[95] + kb_C3bBbH * y[48] - kf_C3bH * y[4] * y[95] +
    kb_C3bH * y[62] - kf_C3bH * y[6] * y[95] + kb_C3bBbH * y[67] - kf_C3bH * y[51] * y[95] + kb_C3bH_ho * y[59] - 
    kf_C3bH * y[53] * y[95] + kb_C3bBbH * y[54] - kf_C3bH * y[69] * y[95] + kb_C3bBbH * y[70] + 
    kcat_C3bH * y[92] * y[59]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]) + 
    kcat_C3bH * y[92] * y[62]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24])) 


    dydt[96] = (-kf_C3bCR1 * y[1] * y[96] + kb_C3bCR1 * y[47] - 
    kf_C3bCR1 * y[4] * y[96] + kb_C3bCR1 * y[64] - 
    kf_iC3bCR1 * y[63] * y[96] + kb_iC3bCR1 * y[65] - 
    kf_C3bCR1 * y[6] * y[96] + kb_C3bBbCR1 * y[68] - 
    kf_C3bCR1 * y[11] * y[96] + kb_C3bCR1 * y[22] - 
    kf_iC3bCR1 * y[23] * y[96] + kb_iC3bCR1 * y[24] - 
    kf_C3bCR1 * y[51] * y[96] + kb_C3bCR1 * y[57] - 
    kf_iC3bCR1 * y[58] * y[96] + kb_iC3bCR1 * y[60] -
    kf_C3bCR1 * y[53] * y[96] + kb_C3bBbCR1 * y[55] - 
    kf_C3bCR1 * y[69] * y[96] + kb_C3bBbCR1 * y[71] +
    kcat_C3bH * y[92] * y[22]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24])+
    kcat_C3bH * y[92] * y[24]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]) +
    kcat_C3bH * y[92] * y[57]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]) + 
    kcat_C3bH * y[92] * y[60]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]) + 
    kcat_C3bH * y[92] * y[64]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]) + 
    kcat_C3bH * y[92] * y[65]/(km_C3bH + y[62] + y[64] + y[65] + y[59] + y[57] + y[60] + y[22] + y[24]) ) 

    dydt[97] = (-kf_C3bBbDAF * y[53] * y[97] + kb_C3bBbDAF * y[56] - 
    kf_C3bBbDAF * y[69] * y[97] + kb_C3bBbDAF * y[72]) 

    dydt[98] = (-kf_C3bBbC3bC5 * y[69] * y[98] + kb_C3bBbC3bC5 * y[73] - 
    kf_C3bBbC3bC5 * y[27] * y[98] + kb_C3bBbC3bC5 * y[28] - 
    kf_C3bBbC3bC5 * y[31] * y[98] + kb_C3bBbC3bC5 * y[33] - 
    kf_C3bBbC3bC5 * y[32] * y[98] + kb_C3bBbC3bC5 * y[34]-
    
    kcat_C3bBbC3b * y[88] * y[28]/(km_C3bBbC3b+ y[98]) - 
    kcat_C3bBbC3b * y[88] * y[33]/(km_C3bBbC3b+ y[98]) - 
    kcat_C3bBbC3b * y[88] * y[34]/(km_C3bBbC3b+ y[98]) - 
    kcat_C3bBbC3b * y[88] * y[73]/(km_C3bBbC3b+ y[98])) 

    dydt[99] = (kcat_C3bBbC3b * y[88] * y[28]/(km_C3bBbC3b+ y[98]) + 
    kcat_C3bBbC3b * y[88] * y[33]/(km_C3bBbC3b+ y[98]) + 
    kcat_C3bBbC3b * y[88] * y[34]/(km_C3bBbC3b+ y[98]) + 
    kcat_C3bBbC3b * y[88] * y[73]/(km_C3bBbC3b+ y[98])) 

    dydt[100] = (+kb_C5b * y[29]
    +kb_C5b * y[74]
    +kb_C5b * y[35]
    +kb_C5b * y[36] + 
    kcat_C3bBbC3b * y[88] * y[28]/(km_C3bBbC3b+ y[98]) + 
    kcat_C3bBbC3b * y[88] * y[33]/(km_C3bBbC3b+ y[98]) + 
    kcat_C3bBbC3b * y[88] * y[34]/(km_C3bBbC3b+ y[98]) + 
    kcat_C3bBbC3b * y[88] * y[73]/(km_C3bBbC3b+ y[98])) 

    dydt[101] = (-kf_C3bBbC3bC5bC6 * y[29] * y[101] + kb_C3bBbC3bC5bC6 * y[30]
    -kf_C3bBbC3bC5bC6 * y[35] * y[101] + kb_C3bBbC3bC5bC6 * y[37]
    -kf_C3bBbC3bC5bC6 * y[36] * y[101] + kb_C3bBbC3bC5bC6 * y[38]
    -kf_C3bBbC3bC5bC6 * y[74] * y[101] + kb_C3bBbC3bC5bC6 * y[75])  

    dydt[102] = (-kf_C5b7 * y[30] * y[102]
    -kf_C5b7 * y[38] * y[102]
    -kf_C5b7 * y[37] * y[102] 
    -kf_C5b7 * y[75] * y[102] + kb_C5b7 * y[39] + kb_C5b7 * y[76])  

    dydt[103] = ( - kf_C5b8 * y[39] * y[103]  - kf_C5b8 * y[76] * y[103]  + kb_C5b8 * y[41] - kf_CnC5b8 * y[85] * y[103] +
    kb_CnC5b8 * y[86]- kf_VnC5b8 * y[82] * y[103] + kb_VnC5b8 * y[83] - kf_C5b8 * y[77] * y[103] -
    kf_C5b8 * y[43] * y[103])  

    dydt[104] = (-kf_C5b9 * y[41] * y[104] + kb_C5b9 * y[42]- kf_CnC5b9 * y[86] * y[104] + kb_CnC5b9 * y[87] -
    kf_VnC5b9 * y[83] * y[104] + kb_VnC5b9 * y[84] - kf_C5b9 * y[44] * y[104] - kf_C5b9 * y[78] * y[104] -
    kf_C5b9 * y[79] * y[104])  

    dydt[105] = ( - kf_VnC5b7 * y[76] * y[105] - kf_VnC5b7 * y[39] * y[105] + kb_VnC5b7 * y[82])  

    dydt[106] = ( - kf_CnC5b7 * y[76] * y[106] - kf_CnC5b7 * y[39] * y[106] + kb_CnC5b7 * y[85])   

    dydt[107] = (- kf_CD59C5b9 * y[79] * y[107] + kb_CD59C5b9 * y[81]) 

end