
function [ParameterSet] = KX_ManualCuration_BOTH(ParameterSet)

    % For TKT2 (same values as in TKT1)
    ParameterSet(47) = ParameterSet(33); % KmF6P
    ParameterSet(48) = ParameterSet(34); % KmGAP
    ParameterSet(49) = ParameterSet(35); % KmE4P
    ParameterSet(50) = ParameterSet(36); % KmXu5P
    ParameterSet(51) = ParameterSet(37); % KmS7P
    ParameterSet(52) = ParameterSet(38); % KmR5P

    ParameterSet(56) = ParameterSet(42); % KiATP
    ParameterSet(57) = ParameterSet(43); % KiCIT
    ParameterSet(58) = ParameterSet(44); % KiGLX
    ParameterSet(59) = ParameterSet(45); % KiRuBP


    % For FBA (same values as in ALD)
    ParameterSet(70) = ParameterSet(26); % KmFBP
    ParameterSet(71) = ParameterSet(27); % KmDHAP
    ParameterSet(72) = ParameterSet(28); % KmGAP
    ParameterSet(73) = ParameterSet(29); % KmSBP
    ParameterSet(74) = ParameterSet(30); % KmE4P

    
    % For SBPase (same values as FBPase)
    ParameterSet(77) = ParameterSet(61); % KmFBP
    ParameterSet(78) = ParameterSet(62); % KmSBP
    ParameterSet(81) = ParameterSet(65); % KiAMP
    ParameterSet(82) = ParameterSet(66); % KiNADPH
    ParameterSet(83) = ParameterSet(67); % KiCIT
    ParameterSet(84) = ParameterSet(68); % KaGAP
    
    
    % For XFPK2 (same values as XFPK1)
    ParameterSet(136) = ParameterSet(128); % KmF6P
    ParameterSet(137) = ParameterSet(129); % KmP_i
    ParameterSet(138) = ParameterSet(130); % KmE4P
    ParameterSet(139) = ParameterSet(131); % KmACETP
    ParameterSet(140) = ParameterSet(132); % KmXu5P
    ParameterSet(141) = ParameterSet(133); % KmGAP


    % For Rubiso (same values as Rubisc)
    ParameterSet(234) = ParameterSet(2); % KmO2
    ParameterSet(235) = ParameterSet(3); % KmCO2_cax
    ParameterSet(236) = ParameterSet(4); % KmRuBP
    ParameterSet(237) = ParameterSet(5); % KaP_i
    ParameterSet(238) = ParameterSet(6); % KiP_i
    ParameterSet(239) = ParameterSet(7); % KiNADPH


end