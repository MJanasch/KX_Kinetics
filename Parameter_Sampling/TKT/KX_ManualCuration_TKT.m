
function [ParameterSet] = KX_ManualCuration_TKT(ParameterSet)

    % For TKT2 (same values as in TKT1)
    ParameterSet(42) = ParameterSet(33); % KmF6P
    ParameterSet(43) = ParameterSet(34); % KmGAP
    ParameterSet(44) = ParameterSet(35); % KmE4P
    ParameterSet(45) = ParameterSet(36); % KmXu5P
    ParameterSet(46) = ParameterSet(37); % KmS7P
    ParameterSet(47) = ParameterSet(38); % KmR5P

    ParameterSet(56) = ParameterSet(42); % KiATP
    ParameterSet(57) = ParameterSet(43); % KiCIT
    ParameterSet(58) = ParameterSet(44); % KiGLX
    ParameterSet(59) = ParameterSet(45); % KiRuBP


    % For FBA (same values as in ALD)
    ParameterSet(66) = ParameterSet(26); % KmFBP
    ParameterSet(67) = ParameterSet(27); % KmDHAP
    ParameterSet(68) = ParameterSet(28); % KmGAP
    ParameterSet(69) = ParameterSet(29); % KmSBP
    ParameterSet(70) = ParameterSet(30); % KmE4P

    
    % For SBPase (same values as FBPase)
    ParameterSet(73) = ParameterSet(61); % KmFBP
    ParameterSet(74) = ParameterSet(62); % KmSBP

    % For XFPK2 (same values as XFPK1)
    ParameterSet(128) = ParameterSet(120); % KmF6P
    ParameterSet(129) = ParameterSet(121); % KmP_i
    ParameterSet(130) = ParameterSet(122); % KmE4P
    ParameterSet(131) = ParameterSet(123); % KmACETP
    ParameterSet(132) = ParameterSet(124); % KmXu5P
    ParameterSet(133) = ParameterSet(125); % KmGAP


    % For Rubiso (same values as Rubisc)
    ParameterSet(226) = ParameterSet(2); % KmO2
    ParameterSet(227) = ParameterSet(3); % KmCO2_cax
    ParameterSet(228) = ParameterSet(4); % KmRuBP
    ParameterSet(229) = ParameterSet(5); % KaP_i
    ParameterSet(230) = ParameterSet(6); % KiP_i
    ParameterSet(231) = ParameterSet(7); % KiNADPH


end