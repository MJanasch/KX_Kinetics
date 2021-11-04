
function [ParameterSet] = KX_ManualCuration_FSBPase(ParameterSet)

    % For TKT2 (same values as in TKT1)
    ParameterSet(42) = ParameterSet(33); % KmF6P
    ParameterSet(43) = ParameterSet(34); % KmGAP
    ParameterSet(44) = ParameterSet(35); % KmE4P
    ParameterSet(45) = ParameterSet(36); % KmXu5P
    ParameterSet(46) = ParameterSet(37); % KmS7P
    ParameterSet(47) = ParameterSet(38); % KmR5P

    % For FBA (same values as in ALD)
    ParameterSet(60) = ParameterSet(26); % KmFBP
    ParameterSet(61) = ParameterSet(27); % KmDHAP
    ParameterSet(62) = ParameterSet(28); % KmGAP
    ParameterSet(63) = ParameterSet(29); % KmSBP
    ParameterSet(64) = ParameterSet(30); % KmE4P

    
    % For SBPase (same values as FBPase)
    ParameterSet(67) = ParameterSet(51); % KmFBP
    ParameterSet(68) = ParameterSet(52); % KmSBP
    ParameterSet(71) = ParameterSet(55); % KiAMP
    ParameterSet(72) = ParameterSet(56); % KiNADPH
    ParameterSet(73) = ParameterSet(57); % KiCIT
    ParameterSet(74) = ParameterSet(58); % KaGAP
    
    
    % For XFPK2 (same values as XFPK1)
    ParameterSet(126) = ParameterSet(118); % KmF6P
    ParameterSet(127) = ParameterSet(119); % KmP_i
    ParameterSet(128) = ParameterSet(120); % KmE4P
    ParameterSet(129) = ParameterSet(121); % KmACETP
    ParameterSet(130) = ParameterSet(122); % KmXu5P
    ParameterSet(131) = ParameterSet(123); % KmGAP


    % For Rubiso (same values as Rubisc)
    ParameterSet(224) = ParameterSet(2); % KmO2
    ParameterSet(225) = ParameterSet(3); % KmCO2_cax
    ParameterSet(226) = ParameterSet(4); % KmRuBP
    ParameterSet(227) = ParameterSet(5); % KaP_i
    ParameterSet(228) = ParameterSet(6); % KiP_i
    ParameterSet(229) = ParameterSet(7); % KiNADPH


end