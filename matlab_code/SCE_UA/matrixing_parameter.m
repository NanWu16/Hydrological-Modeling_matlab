function TotalParam = matrixing_parameter(params)
    index_num   = (1 : 12)';
    OC          = params(1) * ones(12, 1);
    ROC         = params(2) * ones(12, 1);
    
    KEpC        = params(3  : 3  + 11)';
    DeeperC     = params(15 : 15 + 11)';
    AlUpper     = params(27 : 27 + 11)';
    AlLower     = params(39 : 39 + 11)';
    CCg         = params(51 : 51 + 11)';
    CCi         = params(63 : 63 + 11)';
    CCS         = params(75 : 75 + 11)';
    
    LagTime     = params(87) * ones(12, 1);
    cshm        = params(88) * ones(12, 1);
    lagtimehm   = params(89) * ones(12, 1);
    
    MKch        = params(90 : 90 + 11)';
    MKs         = params(102 : 102 + 11)';
    MKi         = params(114 : 114 + 11)';
    MKg         = params(126 : 126 + 11)';
    
    MXch        = params(138) * ones(12, 1);
    MXs         = params(139) * ones(12, 1);
    MXi         = params(140) * ones(12, 1);
    MXg         = params(141) * ones(12, 1);
    
    UADJ        = params(142) * ones(12, 1);
    MBASE       = params(143) * ones(12, 1);
    MFMAX       = params(144) * ones(12, 1);
    MFMIN       = params(145) * ones(12, 1);

    TIPM        = params(146) * ones(12, 1);
    NMF         = params(147) * ones(12, 1);
    PLWHC       = params(148) * ones(12, 1);

    DAYGM       = params(149 : 149 + 11)';
    
    R1 = params(161 : 161 + 11)';
    
    TotalParam = [index_num, OC, ROC, KEpC, DeeperC, AlUpper, AlLower, CCg, CCi, CCS, LagTime, cshm, lagtimehm,...
        MKch, MKs, MKi, MKg, MXch, MXs, MXi, MXg, UADJ, MBASE, MFMAX, MFMIN, TIPM, NMF, PLWHC, DAYGM, R1];
end