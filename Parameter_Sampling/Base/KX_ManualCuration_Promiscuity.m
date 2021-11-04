
function [ParameterSet] = KX_ManualCuration_Promiscuity(ParameterSet)

    % For TKT2 (same values as in TKT1)
    ParameterSet(42) = ParameterSet(33); % KmF6P
    ParameterSet(43) = ParameterSet(34); % KmGAP
    ParameterSet(44) = ParameterSet(35); % KmE4P
    ParameterSet(45) = ParameterSet(36); % KmXu5P
    ParameterSet(46) = ParameterSet(37); % KmS7P
    ParameterSet(47) = ParameterSet(38); % KmR5P

    % For FBA (same values as in ALD)
    ParameterSet(56) = ParameterSet(26); % KmFBP
    ParameterSet(57) = ParameterSet(27); % KmDHAP
    ParameterSet(58) = ParameterSet(28); % KmGAP
    ParameterSet(59) = ParameterSet(29); % KmSBP
    ParameterSet(60) = ParameterSet(30); % KmE4P

    % For SBPase (same values as FBPase)
    ParameterSet(63) = ParameterSet(51); % KmFBP
    ParameterSet(64) = ParameterSet(52); % KmSBP

    % For XFPK2 (same values as XFPK1)
    ParameterSet(118) = ParameterSet(110); % KmF6P
    ParameterSet(119) = ParameterSet(111); % KmP_i
    ParameterSet(120) = ParameterSet(112); % KmE4P
    ParameterSet(121) = ParameterSet(113); % KmACETP
    ParameterSet(122) = ParameterSet(114); % KmXu5P
    ParameterSet(123) = ParameterSet(115); % KmGAP


    % For Rubiso (same values as Rubisc)
    ParameterSet(216) = ParameterSet(2); % KmO2
    ParameterSet(217) = ParameterSet(3); % KmCO2_cax
    ParameterSet(218) = ParameterSet(4); % KmRuBP
    ParameterSet(219) = ParameterSet(5); % KaP_i
    ParameterSet(220) = ParameterSet(6); % KiP_i
    ParameterSet(221) = ParameterSet(7); % KiNADPH


end