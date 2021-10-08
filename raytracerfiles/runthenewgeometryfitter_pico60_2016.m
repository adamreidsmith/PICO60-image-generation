
%% First locate fiducial marks
% This from e-log entry 6097
% http://dbweb5.fnal.gov:8080/ECL/coupp/E/show?e=6097

% 1st dim:  row number
% 2nd dim: column number, column 2 at phi=0, column 4 at phi=pi
% 3rd dim: z, r*phi

fidmarks_zphi = zeros(6,4,2);

fidmarks_zphi(:,:,1) = 105.2 - [ ...
    35.5, 35.6, 35.5, 35.5 ; ...
    47, 47, 47.3, 47 ; ...
    55, 55, 55, 55 ; ...
    64, 64.2, 64, 64 ; ...
    75, 75, 75, 75 ; ...
    86, 86, 86, 86];
    
fidmarks_zphi(:,4,2) = 47.2;
fidmarks_zphi(:, 1, 2) = ...
    -[32.25 ; 32.25 ; 32.5 ; 32.5 ; 32.8 ; 33];
fidmarks_zphi(:, 3, 2) = ...
    [32.75 ; 33.75 ; 34 ; 34.5 ; 34.7 ; 35.5];

%% Now start some izphi's and pixel lists -- first fiducial marks
% this is cam#, ipix, jpix, markrow, markcolumn, marktype
% pixels here are (1,1) corresponds to upper-RIGHT of rotated image, ipix
% is vertical coordinate in rotated image, jpix is horizontal coordinate in
% rotated image (i.e. I'm grabbing pixels off un-rotated images) )
% ipix in 1:1700, jpix in 1:1088
% marktype can be 0 (normal), 1 (z-only), or 2 (phi-only)
fidmarks_pix = [ ...
    1, 577, 425, 1, 3, 0 ; ... 1:3
    1, 577, 417, 1, 3, 1 ; ... 4
    1, 572, 425, 1, 3, 2 ; ... 5:6
    1, 582, 425, 1, 3, 2 ; ... 7:8
    1, 846, 420, 2, 3, 0 ; ... 9:11
    1, 846, 415, 2, 3, 1 ; ... 12
    1, 843, 420, 2, 3, 2 ; ... 13:14
    1, 849, 420, 2, 3, 2 ; ... 15:16
    1, 1033, 418, 3, 3, 0 ; ... 17:19
    1, 1033, 414, 3, 3, 1 ; ... 20
    1, 1029, 419, 3, 3, 2 ; ... 21:22
    1, 1036, 418, 3, 3, 2 ; ... 23:24
    1, 1236, 416, 4, 3, 0 ; ... 25:27
    1, 1236, 411, 4, 3, 1 ; ... 28
    1, 1232, 417, 4, 3, 2 ; ... 29:30
    1, 1240, 416, 4, 3, 2 ; ... 31:32
    1, 1463, 419, 5, 3, 0 ; ... 33:35
    1, 1463, 414, 5, 3, 1 ; ... 36
    1, 1458, 420, 5, 3, 2 ; ... 37:38
    1, 1662, 427, 6, 3, 0 ; ... 39:41
    1, 570, 714, 1, 4, 0 ; ... 42
    1, 846, 711, 2, 4, 0 ; ... 43
    1, 1038, 702, 3, 4, 0 ; ... 44
    1, 1246, 692, 4, 4, 0 ; ... 45
    1, 1477, 672, 5, 4, 0 ; ... 46
    1, 508, 924, 1, 1, 0 ; ... 47:49
    1, 511, 921, 1, 1, 1 ; ... 50
    1, 503, 924, 1, 1, 2 ; ... 51:52
    1, 513, 924, 1, 1, 2 ; ... 53:54
    1, 872, 925, 2, 1, 0 ; ... 55:57
    1, 872, 922, 2, 1, 1 ; ... 58
    1, 867, 925, 2, 1, 2 ; ... 59:60
    1, 877, 925, 2, 1, 2 ; ... 61:62
    1, 1122, 913, 3, 1, 0 ; ... 63:65
    1, 1119, 908, 3, 1, 1 ; ... 66
    1, 1115, 913, 3, 1, 2 ; ... 67:68
    1, 1127, 913, 3, 1, 2 ; ... 69:70
    1, 1385, 882, 4, 1, 0 ; ... 71:73
    1, 1382, 880, 4, 1, 1 ; ... 74
    1, 1381, 883, 4, 1, 2 ; ... 75:76
    1, 1390, 882, 4, 1, 2 ; ... 77:78
    1, 318, 134, 1, 2, 0 ; ... 71:73 +8
    1, 871, 110, 2, 2, 0 ; ... 74:76 +8
    1, 1259, 106, 3, 2, 0 ; ... 77:79 +8
    2, 143, 465, 1, 3, 0 ; ... 80:82 +8
    2, 143, 461, 1, 3, 1 ; ... 83 +8
    2, 146, 465, 1, 3, 2 ; ... 84:85 +8
    2, 333, 456, 2, 3, 0 ; ... 86:88 +8
    2, 481, 451, 3, 3, 0 ; ... 89:91 +8
    2, 659, 445, 4, 3, 0 ; ... 92:94 +8
    2, 889, 444, 5, 3, 0 ; ... 95:97 +8
    2, 1124, 450, 6, 3, 0 ; ... 98:100 +8
    2, 1124, 446, 6, 3, 1 ; ... 101 +8
    2, 125, 682, 1, 4, 0 ; ... 102 +8
    2, 317, 691, 2, 4, 0 ; ... 103 +8
    2, 467, 696, 3, 4, 0 ; ... 104 +8
    2, 650, 700, 4, 4, 0 ; ... 105 +8
    2, 885, 700, 5, 4, 0 ; ... 106 +8
    2, 1125, 694, 6, 4, 0 ; ... 107 +8
    2, 179, 849, 2, 1, 0 ; ... 108:110 +8
    2, 362, 870, 3, 1, 0 ; ... 111:113 +8
    2, 367, 867, 3, 1, 1 ; ... 114 +8
    2, 355, 870, 3, 1, 2 ; ... 115:116 +8
    2, 596, 887, 4, 1, 0 ; ... 117:119 +8
    2, 904, 894, 5, 1, 0 ; ... 120:122 +8
    2, 898, 894, 5, 1, 2 ; ... 123:124 +8
    2, 1208, 882, 6, 1, 0 ; ... 125:127 +8
    2, 440, 171, 4, 2, 0 ; ... 128:130 +8
    2, 923, 154, 5, 2, 0 ; ... 131:133 +8
    2, 1391, 162, 6, 2, 0 ; ... 134:136 +8
    3, 611, 707, 1, 1, 0 ; ... 137:139 +8
    3, 873, 710, 2, 1, 0 ; ... 140:142 +8
    3, 873, 704, 2, 1, 1 ; ... 143 +8
    3, 867, 710, 2, 1, 2 ; ... 144:145 +8
    3, 877, 710, 2, 1, 2 ; ... 146:147 +8
    3, 1054, 709, 3, 1, 0 ; ... 148:150 +8
    3, 1054, 696, 3, 1, 1 ; ... 151 +8
    3, 1049, 709, 3, 1, 2 ; ... 152:153 +8
    3, 1058, 709, 3, 1, 2 ; ... 154:155 +8
    3, 1249, 703, 4, 1, 0 ; ... 156:158 +8
    3, 1249, 698, 4, 1, 2 ; ... 159:160 +8
    3, 1467, 692, 5, 1, 0 ; ... 161:163 +8
    3, 1655, 678, 6, 1, 0 ; ... 164:166 +8
    3, 355, 873, 1, 2, 0 ; ... 167:169 +8
    3, 904, 887, 2, 2, 0 ; ... 170:172 +8
    3, 1289, 878, 3, 2, 0 ; ... 173:175 +8
    3, 550, 215, 1, 3, 0 ; ... 176:178 +8
    3, 874, 210, 2, 3, 0 ; ... 179:181 +8
    3, 1100, 218, 3, 3, 0 ; ... 182:184 +8
    3, 1337, 235, 4, 3, 0 ; ... 185:187 +8
    3, 1591, 266, 5, 3, 0 ; ... 188:190 +8
    3, 604, 432, 1, 4, 0 ; ... 191 +8
    3, 862, 432, 2, 4, 0 ; ... 192 +8
    3, 1042, 432, 3, 4, 0 ; ... 193 +8
    3, 1236, 437, 4, 4, 0 ; ... 194 +8
    3, 1453, 444, 5, 4, 0 ; ... 195 +8
    4, 128, 682, 1, 1, 0 ; ... 196:198 +8
    4, 320, 686, 2, 1, 0 ; ... 199:201 +8
    4, 321, 682, 2, 1, 1 ; ... 202 +8
    4, 317, 686, 2, 1, 2 ; ... 203:204 +8
    4, 469, 688, 3, 1, 0 ; ... 205:207 +8
    4, 470, 681, 3, 1, 1 ; ... 208 +8
    4, 464, 688, 3, 1, 2 ; ... 209:210 +8
    4, 650, 685, 4, 1, 0 ; ... 211:213 +8
    4, 651, 680, 4, 1, 1 ; ... 214 +8
    4, 646, 685, 4, 1, 2 ; ... 215:216 +8
    4, 881, 676, 5, 1, 0 ; ... 217:219 +8
    4, 881, 671, 5, 1, 1 ; ... 220 +8
    4, 876, 676, 5, 1, 2 ; ... 221:222 +8
    4, 884, 676, 5, 1, 2 ; ... 223:224 +8
    4, 1112, 662, 6, 1, 0 ; ... 225:227 +8
    4, 76, 834, 3, 2, 0 ; ... 228:230 +8
    4, 430, 839, 4, 2, 0 ; ... 231:233 +8
    4, 921, 824, 5, 2, 0 ; ... 234:236 +8
    4, 1394, 794, 6, 2, 0 ; ... 237:239 +8
    4, 186, 285, 2, 3, 0 ; ... 240:242 +8
    4, 361, 261, 3, 3, 0 ; ... 243:245 +8
    4, 577, 237, 4, 3, 0 ; ... 246:248 +8
    4, 861, 222, 5, 3, 0 ; ... 249:251 +8
    4, 1144, 224, 6, 3, 0 ; ... 252:254 +8
    4, 122, 470, 1, 4, 0 ; ... 255 +8
    4, 309, 454, 2, 4, 0 ; ... 256 +8
    4, 456, 444, 3, 4, 0 ; ... 257 +8
    4, 633, 433, 4, 4, 0 ; ... 258 +8
    4, 859, 423, 5, 4, 0 ; ... 259 +8
    4, 1091, 413, 6, 4, 0 ; ... 260 +8
    ];
    
fid_izphi = zeros(size(fidmarks_pix, 1), 3);
fid_izphi(:, 1) = 4;
fid_izphi(fidmarks_pix(:,5)==2, 1) = 1;
fid_izphi(:, 2) = diag( ...
    fidmarks_zphi(fidmarks_pix(:, 4), fidmarks_pix(:, 5), 1));
fid_izphi(:, 3) = diag( ...
    fidmarks_zphi(fidmarks_pix(:, 4), fidmarks_pix(:, 5), 2));
fid_izphi(fidmarks_pix(:,6)==1, 3) = NaN;
fid_izphi(fidmarks_pix(:,6)==2, 2) = NaN;
fid_izphi(fidmarks_pix(:,5)==4, 3) = NaN;
fid_pix = fidmarks_pix(:, 1:3);
fid_err = zeros(size(fidmarks_pix,1),1) + 0.3;

%% Now some other types of points    
%%% apex of the jar
apex_pix = [ ...
    2, 1425, 557, ; ... 261:263
    4, 1385, 524, ; ... 264:266
    ];
apex_izphi = repmat([4, -inf, nan], size(apex_pix,1), 1);
apex_err = zeros(size(apex_pix,1), 1) + 0.5;

%%% front of the target-buffer interface
frontsurface_pix = [ ...
    1, 256, 944 ; ...
    1, 202, 766 ; ...
    1, 180, 557 ; ...
    1, 202, 282 ; ...
    1, 249, 129 ; ...
    3, 242, 561 ; ...
    3, 282, 263 ; ...
    ];
frontsurface_izphi = repmat([2, inf, nan], size(frontsurface_pix,1), 1);
frontsurface_err = zeros(size(frontsurface_pix,1), 1) + 0.5;

%%% back of the target-buffer interface
backsurface_pix = [ ...
    1, 448, 925 ; ...
    1, 534, 752 ; ...
    1, 550, 588 ; ...
    1, 529, 349 ; ...
    1, 460, 191 ; ...
    2, 99, 418 ; ...
    2, 99, 691 ; ...
    2, 130, 590 ; ...
    3, 532, 255 ; ...
    3, 576, 425 ; ...
    3, 582, 571 ; ...
    3, 566, 736 ; ...
    4, 93, 695 ; ...
    4, 119, 572 ; ...
    4, 84, 435 ; ...
    ];
backsurface_izphi = repmat([3, inf, nan], size(backsurface_pix,1), 1);
backsurface_err = zeros(size(backsurface_pix,1), 1) + 0.2;

%%% tangent to outside of jar
odtan_pix = [ ...
    1, 256, 13 ; ...
    1, 1206, 1075 ; ...
    1, 1448, 1057 ; ...
    1, 1092, 1082 ; ...
    1, 1311, 1068 ; ...
    1, 325, 7 ; ...
    2, 1615, 558 ; ...
    2, 1129, 53 ; ...
    2, 384, 69 ; ...
    2, 375, 1031 ; ...
    2, 1105, 1042 ; ...
    2, 328, 74 ; ...
    2, 601, 58 ; ...
    2, 811, 52 ; ...
    2, 1019, 51 ; ...
    2, 1326, 62 ; ...
    2, 1479, 133 ; ...
    2, 1595, 370 ; ...
    2, 1594, 736 ; ...
    2, 1539, 883 ; ...
    2, 1439, 999 ; ...
    2, 1255, 1036 ; ...
    2, 927, 1042 ; ...
    2, 1556, 250 ; ...
    2, 728, 1042 ; ...
    2, 532, 1036 ; ...
    3, 310, 47 ; ...
    3, 707, 36 ; ...
    3, 1065, 40 ; ...
    3, 1407, 55 ; ...
    3, 531, 39 ; ...
    3, 881, 38 ; ...  
    3, 1198, 45 ; ...
    3, 403, 1069 ; ...
    3, 892, 1085 ; ...
    3, 1400, 1071 ; ...
    3, 1270, 1078 ; ...
    3, 1045, 1085 ; ...
    3, 723, 1083 ; ...
    3, 513, 1075 ; ...
    3, 247, 1058 ; ...
    ...4, 312, 102 ; ...
    4, 563, 75 ; ...
    4, 800, 63 ; ...
    4, 1010, 59 ; ...
    4, 1423, 126 ; ...
    4, 1510, 264 ; ...
    4, 1553, 427 ; ...
    4, 1552, 686 ; ...
    4, 1516, 819 ; ...
    4, 1453, 912 ; ...
    4, 1259, 1000 ; ...
    4, 994, 1014 ; ...
    4, 800, 1023 ; ...
    4, 596, 1025 ; ...
    4, 1563, 554 ; ...
    4, 1218, 59 ; ...
    4, 420, 87 ; ...
    4, 1128, 1008 ; ...
    4, 399, 1025 ; ...
    ];
odtan_izphi = repmat([-1, nan, nan], size(odtan_pix,1), 1);
odtan_err = zeros(size(odtan_pix,1), 1) + 4;

%%% tangent to inside of jar
idtan_pix = [ ...
    1, 1217, 1065 ; ...
    2, 404, 1019 ; ...
    2, 1595, 592 ; ...
    2, 1130, 60 ; ...
    2, 380, 79 ; ...
    2, 1156, 1027 ; ...
    3, 692, 46 ; ...
    3, 1221, 54 ; ...
    3, 729, 1070 ; ...
    3, 1325, 1062 ; ...
    4, 1095, 68 ; ...
    4, 373, 101 ; ...
    4, 1136, 1000 ; ...
    4, 397, 1016 ; ...
    4, 1543, 534 ; ...
    ];
idtan_izphi = repmat([-3, nan, nan], size(idtan_pix,1), 1);
idtan_err = zeros(size(idtan_pix,1), 1) + 4;

%%% edge of target volume
crit_pix = [ ...
    1, 1424, 133 ; ...
    1, 1020, 109 ; ...
    1, 493, 129 ; ...
    1, 1044, 979 ; ...
    1, 575, 987 ; ...
    2, 347, 918 ; ...
    2, 868, 163 ; ...
    2, 1524, 582 ; ...
    2, 1234, 927 ; ...
    3, 726, 130 ; ...
    3, 546, 950 ; ...
    3, 898, 963 ; ...
    3, 1308, 942 ; ...
    3, 1123, 141 ; ...
    4, 1181, 887 ; ...
    4, 544, 911 ; ...
    4, 866, 913 ; ...
    4, 1480, 534 ; ...
    4, 1156, 161 ; ...
    4, 234, 215 ; ...
    4, 728, 166 ; ...
    ];
crit_izphi = repmat([-4, nan, nan], size(crit_pix,1), 1);
crit_err = zeros(size(crit_pix,1), 1) + 0.1;

%% now make our initial params

cam = struct('x',[],'y',[],'z',[],'pitch',[],'yaw',[],'roll',[], ...
    'f',[],'f1',[],'f2',[],'i0',[],'j0',[],'i_pitch',[],'j_pitch',[],'theta',[],'phi',[],'bf',[],'lens_type',[], ...
    'win_d',[],'win_t',[],'win_phi',[],'win_pitch',[],'win_yaw',[]);

starting_params = 'secondtry';

freshfit = 1;

if ~exist('paramstruct','var') || ~isstruct(paramstruct) || isempty(paramstruct)
    paramstruct = struct();
end

paramsets = {'lowerbound','upperbound',starting_params};
for np = 1:length(paramsets)
    switch paramsets{np}
        case 'nominal'
            cam(1).x = 0;
            cam(1).y = -1.6;
            cam(1).z = 3.5 + 13*2.54;
            cam(1).pitch = 0;
            cam(1).yaw = -pi/6;
            cam(1).roll = pi/2;
            cam(1).f = .5;
            cam(1).f1 = 0;
            cam(1).f2 = 0;
            cam(1).i0 = .5*2049 - 170;
            cam(1).j0 = .5*1089;
            cam(1).i_pitch = .00055;
            cam(1).j_pitch = .00055;
            cam(1).theta = 0;
            cam(1).phi = 0;
            cam(1).bf = 10;
            cam(1).lens_type = 1;
            cam(1).win_d = 15*2.54;
            cam(1).win_t = 2*2.54;
            cam(1).win_phi = -pi/6;
            cam(1).win_pitch = 0;
            cam(1).win_yaw = -pi/6;

            cam(2).x = 0;
            cam(2).y = -1.6;
            cam(2).z = 3.5 + 2*2.54;
            cam(2).pitch = 0;
            cam(2).yaw = -pi/6;
            cam(2).roll = pi/2;
            cam(2).f = .5;
            cam(2).f1 = 0;
            cam(2).f2 = 0;
            cam(2).i0 = .5*2049 - 170;
            cam(2).j0 = .5*1089;
            cam(2).i_pitch = .00055;
            cam(2).j_pitch = .00055;
            cam(2).theta = 0;
            cam(2).phi = 0;
            cam(2).bf = 10;
            cam(2).lens_type = 1;
            cam(2).win_d = 15*2.54;
            cam(2).win_t = 2*2.54;
            cam(2).win_phi = -pi/6;
            cam(2).win_pitch = 0;
            cam(2).win_yaw = -pi/6;

            cam(3).x = 0;
            cam(3).y = -1.6;
            cam(3).z = 3.5 + 13*2.54;
            cam(3).pitch = 0;
            cam(3).yaw = pi/6;
            cam(3).roll = pi/2;
            cam(3).f = .5;
            cam(3).f1 = 0;
            cam(3).f2 = 0;
            cam(3).i0 = .5*2049 - 170;
            cam(3).j0 = .5*1089;
            cam(3).i_pitch = .00055;
            cam(3).j_pitch = .00055;
            cam(3).theta = 0;
            cam(3).phi = 0;
            cam(3).bf = 10;
            cam(3).lens_type = 1;
            cam(3).win_d = 15*2.54;
            cam(3).win_t = 2*2.54;
            cam(3).win_phi = pi/6;
            cam(3).win_pitch = 0;
            cam(3).win_yaw = pi/6;

            cam(4).x = 0;
            cam(4).y = -1.6;
            cam(4).z = 3.5 + 2*2.54;
            cam(4).pitch = 0;
            cam(4).yaw = pi/6;
            cam(4).roll = pi/2;
            cam(4).f = .5;
            cam(4).f1 = 0;
            cam(4).f2 = 0;
            cam(4).i0 = .5*2049 - 170;
            cam(4).j0 = .5*1089;
            cam(4).i_pitch = .00055;
            cam(4).j_pitch = .00055;
            cam(4).theta = 0;
            cam(4).phi = 0;
            cam(4).bf = 10;
            cam(4).lens_type = 1;
            cam(4).win_d = 15*2.54;
            cam(4).win_t = 2*2.54;
            cam(4).win_phi = pi/6;
            cam(4).win_pitch = 0;
            cam(4).win_yaw = pi/6;
            
            jar_phi = 4*pi/180;
            liquidlevel = 53;

            jar_wall = .5;
            jar_OD = 30;
            jar_axwall = .5;
            jar_Oaxrad = 30;
            jar_knucklewall = .5;
            jar_knuckleOrad = 5;
            jar_xpitch = 0;
            jar_ypitch = 0;

            n_air = 1;
            n_window = 1.458;
            n_hydraulic = 1.434;
            n_jar = 1.458;
            n_target = 1.197;
            
            z_offset = 7.0871;

        case 'upperbound'
            cam(1).x = 5;
            cam(1).y = -.1;
            cam(1).z = 3.5 + 13*2.54 + 5;
            cam(1).pitch = pi/12;
            cam(1).yaw = -pi/6 + pi/12;
            cam(1).roll = pi/2 + pi/12;
            cam(1).f = .8;
            cam(1).f1 = 1;
            cam(1).f2 = 1;
            cam(1).i0 = .5*2049 - 170;
            cam(1).j0 = .5*1089;
            cam(1).i_pitch = .00055;
            cam(1).j_pitch = .00055;
            cam(1).theta = 0;
            cam(1).phi = 0;
            cam(1).bf = 10;
            cam(1).lens_type = 1;
            cam(1).win_d = 15*2.54;
            cam(1).win_t = 2*2.54;
            cam(1).win_phi = -pi/6;
            cam(1).win_pitch = 0;
            cam(1).win_yaw = -pi/6;

            cam(2).x = 5;
            cam(2).y = -.1;
            cam(2).z = 3.5 + 2*2.54 + 5;
            cam(2).pitch = pi/12;
            cam(2).yaw = -pi/6 + pi/12;
            cam(2).roll = pi/2 + pi/12;
            cam(2).f = .8;
            cam(2).f1 = 1;
            cam(2).f2 = 1;
            cam(2).i0 = .5*2049 - 170;
            cam(2).j0 = .5*1089;
            cam(2).i_pitch = .00055;
            cam(2).j_pitch = .00055;
            cam(2).theta = 0;
            cam(2).phi = 0;
            cam(2).bf = 10;
            cam(2).lens_type = 1;
            cam(2).win_d = 15*2.54;
            cam(2).win_t = 2*2.54;
            cam(2).win_phi = -pi/6;
            cam(2).win_pitch = 0;
            cam(2).win_yaw = -pi/6;

            cam(3).x = 5;
            cam(3).y = -.1;
            cam(3).z = 3.5 + 13*2.54 + 5;
            cam(3).pitch = pi/12;
            cam(3).yaw = pi/6 + pi/12;
            cam(3).roll = pi/2 + pi/12;
            cam(3).f = .8;
            cam(3).f1 = 1;
            cam(3).f2 = 1;
            cam(3).i0 = .5*2049 - 170;
            cam(3).j0 = .5*1089;
            cam(3).i_pitch = .00055;
            cam(3).j_pitch = .00055;
            cam(3).theta = 0;
            cam(3).phi = 0;
            cam(3).bf = 10;
            cam(3).lens_type = 1;
            cam(3).win_d = 15*2.54;
            cam(3).win_t = 2*2.54;
            cam(3).win_phi = pi/6;
            cam(3).win_pitch = 0;
            cam(3).win_yaw = pi/6;

            cam(4).x = 5;
            cam(4).y = -.1;
            cam(4).z = 3.5 + 2*2.54 + 5;
            cam(4).pitch = pi/12;
            cam(4).yaw = pi/6 + pi/12;
            cam(4).roll = pi/2 + pi/12;
            cam(4).f = .8;
            cam(4).f1 = 1;
            cam(4).f2 = 1;
            cam(4).i0 = .5*2049 - 170;
            cam(4).j0 = .5*1089;
            cam(4).i_pitch = .00055;
            cam(4).j_pitch = .00055;
            cam(4).theta = 0;
            cam(4).phi = 0;
            cam(4).bf = 10;
            cam(4).lens_type = 1;
            cam(4).win_d = 15*2.54;
            cam(4).win_t = 2*2.54;
            cam(4).win_phi = pi/6;
            cam(4).win_pitch = 0;
            cam(4).win_yaw = pi/6;
            
            jar_phi = pi/6;
            liquidlevel = 53 + 5;

            jar_wall = .5;
            jar_OD = 30;
            jar_axwall = .5;
            jar_Oaxrad = 30;
            jar_knucklewall = .5;
            jar_knuckleOrad = 5;
            jar_xpitch = 0;
            jar_ypitch = 0;

            n_air = 1;
            n_window = 1.458;
            n_hydraulic = 1.434 + .05;
            n_jar = 1.458;
            n_target = 1.197 + .05;
            
            z_offset = 7.0871 + 20;

        case 'lowerbound'
            cam(1).x = -5;
            cam(1).y = -10;
            cam(1).z = 7.9 + 13*2.54 - 5;
            cam(1).pitch = -pi/12;
            cam(1).yaw = -pi/6 - pi/12;
            cam(1).roll = pi/2 - pi/12;
            cam(1).f = .3;
            cam(1).f1 = -1;
            cam(1).f2 = -1;
            cam(1).i0 = .5*2049 - 170;
            cam(1).j0 = .5*1089;
            cam(1).i_pitch = .00055;
            cam(1).j_pitch = .00055;
            cam(1).theta = 0;
            cam(1).phi = 0;
            cam(1).bf = 10;
            cam(1).lens_type = 1;
            cam(1).win_d = 15*2.54;
            cam(1).win_t = 2*2.54;
            cam(1).win_phi = -pi/6;
            cam(1).win_pitch = 0;
            cam(1).win_yaw = -pi/6;

            cam(2).x = -5;
            cam(2).y = -10;
            cam(2).z = 7.9 + 2*2.54 - 5;
            cam(2).pitch = -pi/12;
            cam(2).yaw = -pi/6 - pi/12;
            cam(2).roll = pi/2 - pi/12;
            cam(2).f = .3;
            cam(2).f1 = -1;
            cam(2).f2 = -1;
            cam(2).i0 = .5*2049 - 170;
            cam(2).j0 = .5*1089;
            cam(2).i_pitch = .00055;
            cam(2).j_pitch = .00055;
            cam(2).theta = 0;
            cam(2).phi = 0;
            cam(2).bf = 10;
            cam(2).lens_type = 1;
            cam(2).win_d = 15*2.54;
            cam(2).win_t = 2*2.54;
            cam(2).win_phi = -pi/6;
            cam(2).win_pitch = 0;
            cam(2).win_yaw = -pi/6;

            cam(3).x = -5;
            cam(3).y = -10;
            cam(3).z = 7.9 + 13*2.54 - 5;
            cam(3).pitch = -pi/12;
            cam(3).yaw = pi/6 - pi/12;
            cam(3).roll = pi/2 - pi/12;
            cam(3).f = .3;
            cam(3).f1 = -1;
            cam(3).f2 = -1;
            cam(3).i0 = .5*2049 - 170;
            cam(3).j0 = .5*1089;
            cam(3).i_pitch = .00055;
            cam(3).j_pitch = .00055;
            cam(3).theta = 0;
            cam(3).phi = 0;
            cam(3).bf = 10;
            cam(3).lens_type = 1;
            cam(3).win_d = 15*2.54;
            cam(3).win_t = 2*2.54;
            cam(3).win_phi = pi/6;
            cam(3).win_pitch = 0;
            cam(3).win_yaw = pi/6;

            cam(4).x = -5;
            cam(4).y = -10;
            cam(4).z = 7.9 + 2*2.54 - 5;
            cam(4).pitch = -pi/12;
            cam(4).yaw = pi/6 - pi/12;
            cam(4).roll = pi/2 - pi/12;
            cam(4).f = .3;
            cam(4).f1 = -1;
            cam(4).f2 = -1;
            cam(4).i0 = .5*2049 - 170;
            cam(4).j0 = .5*1089;
            cam(4).i_pitch = .00055;
            cam(4).j_pitch = .00055;
            cam(4).theta = 0;
            cam(4).phi = 0;
            cam(4).bf = 10;
            cam(4).lens_type = 1;
            cam(4).win_d = 15*2.54;
            cam(4).win_t = 2*2.54;
            cam(4).win_phi = pi/6;
            cam(4).win_pitch = 0;
            cam(4).win_yaw = pi/6;
            
            jar_phi = -pi/6;
            liquidlevel = 53 - 5;

            jar_wall = .5;
            jar_OD = 30;
            jar_axwall = .5;
            jar_Oaxrad = 30;
            jar_knucklewall = .5;
            jar_knuckleOrad = 5;
            jar_xpitch = 0;
            jar_ypitch = 0;

            n_air = 1;
            n_window = 1.458;
            n_hydraulic = 1.434 - .05;
            n_jar = 1.458;
            n_target = 1.197 - .05;
            
            z_offset = 7.0871 - 5;
            
        case 'smartstart'
            cam(1).x = 0;
            cam(1).y = -2.5;
            cam(1).z = 3.5 + 13*2.54;
            cam(1).pitch = 0;
            cam(1).yaw = -pi/6;
            cam(1).roll = pi/2;
            cam(1).f = .5;
            cam(1).f1 = 0;
            cam(1).f2 = 0;
            cam(1).i0 = .5*2049 - 170;
            cam(1).j0 = .5*1089;
            cam(1).i_pitch = .00055;
            cam(1).j_pitch = .00055;
            cam(1).theta = 0;
            cam(1).phi = 0;
            cam(1).bf = 10;
            cam(1).lens_type = 1;
            cam(1).win_d = 15*2.54;
            cam(1).win_t = 2*2.54;
            cam(1).win_phi = -pi/6;
            cam(1).win_pitch = 0;
            cam(1).win_yaw = -pi/6;

            cam(2).x = 0;
            cam(2).y = -2.5;
            cam(2).z = 3.5 + 2*2.54;
            cam(2).pitch = 0;
            cam(2).yaw = -pi/6;
            cam(2).roll = pi/2;
            cam(2).f = .5;
            cam(2).f1 = 0;
            cam(2).f2 = 0;
            cam(2).i0 = .5*2049 - 170;
            cam(2).j0 = .5*1089;
            cam(2).i_pitch = .00055;
            cam(2).j_pitch = .00055;
            cam(2).theta = 0;
            cam(2).phi = 0;
            cam(2).bf = 10;
            cam(2).lens_type = 1;
            cam(2).win_d = 15*2.54;
            cam(2).win_t = 2*2.54;
            cam(2).win_phi = -pi/6;
            cam(2).win_pitch = 0;
            cam(2).win_yaw = -pi/6;

            cam(3).x = 0;
            cam(3).y = -2.5;
            cam(3).z = 3.5 + 13*2.54;
            cam(3).pitch = 0;
            cam(3).yaw = pi/6;
            cam(3).roll = pi/2;
            cam(3).f = .5;
            cam(3).f1 = 0;
            cam(3).f2 = 0;
            cam(3).i0 = .5*2049 - 170;
            cam(3).j0 = .5*1089;
            cam(3).i_pitch = .00055;
            cam(3).j_pitch = .00055;
            cam(3).theta = 0;
            cam(3).phi = 0;
            cam(3).bf = 10;
            cam(3).lens_type = 1;
            cam(3).win_d = 15*2.54;
            cam(3).win_t = 2*2.54;
            cam(3).win_phi = pi/6;
            cam(3).win_pitch = 0;
            cam(3).win_yaw = pi/6;

            cam(4).x = 0;
            cam(4).y = -2.5;
            cam(4).z = 3.5 + 2*2.54;
            cam(4).pitch = 0;
            cam(4).yaw = pi/6;
            cam(4).roll = pi/2;
            cam(4).f = .5;
            cam(4).f1 = 0;
            cam(4).f2 = 0;
            cam(4).i0 = .5*2049 - 170;
            cam(4).j0 = .5*1089;
            cam(4).i_pitch = .00055;
            cam(4).j_pitch = .00055;
            cam(4).theta = 0;
            cam(4).phi = 0;
            cam(4).bf = 10;
            cam(4).lens_type = 1;
            cam(4).win_d = 15*2.54;
            cam(4).win_t = 2*2.54;
            cam(4).win_phi = pi/6;
            cam(4).win_pitch = 0;
            cam(4).win_yaw = pi/6;
            
            jar_phi = 5.5*pi/180;
            liquidlevel = 53;

            jar_wall = .5;
            jar_OD = 30;
            jar_axwall = .5;
            jar_Oaxrad = 30;
            jar_knucklewall = .5;
            jar_knuckleOrad = 5;
            jar_xpitch = 0;
            jar_ypitch = 0;

            n_air = 1;
            n_window = 1.45846;
            n_hydraulic = 1.430;
            n_jar = 1.45846;
            n_target = 1.188528;
            
            z_offset = 7.0871;

        case 'firsttry'
            cam(1).x = 0.7666;%0;
            cam(1).y = -2.8509;%-2.5;
            cam(1).z = 37.0873;%7.9 + 13*2.54;
            cam(1).pitch = 0.0293;%0;
            cam(1).yaw = -0.4791;%-pi/6;
            cam(1).roll = 1.5970;%pi/2;
            cam(1).f = 0.6019;%.5;
            cam(1).f1 = -0.0945;%-0.02;
            cam(1).f2 = 0;
            cam(1).i0 = .5*2049 - 170;
            cam(1).j0 = .5*1089;
            cam(1).i_pitch = .00055;
            cam(1).j_pitch = .00055;
            cam(1).theta = 0;
            cam(1).phi = 0;
            cam(1).bf = 10;
            cam(1).lens_type = 1;
            cam(1).win_d = 15*2.54;
            cam(1).win_t = 2*2.54;
            cam(1).win_phi = -pi/6;
            cam(1).win_pitch = 0;
            cam(1).win_yaw = -pi/6;

            cam(2).x = 1.1768;%0;
            cam(2).y = -2.7983;%-2.5;
            cam(2).z = 13.1732;%3.5 + 2*2.54;
            cam(2).pitch = 0.0237;%0;
            cam(2).yaw = -0.4838;%-pi/6;
            cam(2).roll = 1.573;%pi/2;
            cam(2).f = 0.5149;%.5;
            cam(2).f1 = -0.0108;%0;
            cam(2).f2 = 0;
            cam(2).i0 = .5*2049 - 170;
            cam(2).j0 = .5*1089;
            cam(2).i_pitch = .00055;
            cam(2).j_pitch = .00055;
            cam(2).theta = 0;
            cam(2).phi = 0;
            cam(2).bf = 10;
            cam(2).lens_type = 1;
            cam(2).win_d = 15*2.54;
            cam(2).win_t = 2*2.54;
            cam(2).win_phi = -pi/6;
            cam(2).win_pitch = 0;
            cam(2).win_yaw = -pi/6;

            cam(3).x = 0.4623;%0;
            cam(3).y = -3.2301;%-2.5;
            cam(3).z = 40.3593;%3.5 + 13*2.54;
            cam(3).pitch = 0.0020;%0;
            cam(3).yaw = 0.5247;%pi/6;
            cam(3).roll = 1.5643;%pi/2;
            cam(3).f = 0.5560;%.5;
            cam(3).f1 = -0.0215;%0;
            cam(3).f2 = 0;
            cam(3).i0 = .5*2049 - 170;
            cam(3).j0 = .5*1089;
            cam(3).i_pitch = .00055;
            cam(3).j_pitch = .00055;
            cam(3).theta = 0;
            cam(3).phi = 0;
            cam(3).bf = 10;
            cam(3).lens_type = 1;
            cam(3).win_d = 15*2.54;
            cam(3).win_t = 2*2.54;
            cam(3).win_phi = pi/6;
            cam(3).win_pitch = 0;
            cam(3).win_yaw = pi/6;

            cam(4).x = 0.3921;%0;
            cam(4).y = -1.1974;%-2.5;
            cam(4).z = 12.7449;%3.5 + 2*2.54;
            cam(4).pitch = -0.0021;%0;
            cam(4).yaw = 0.5432;%pi/6;
            cam(4).roll = 1.6034;%pi/2;
            cam(4).f = 0.467;%.5;
            cam(4).f1 = -0.0201;%0;
            cam(4).f2 = 0;
            cam(4).i0 = .5*2049 - 170;
            cam(4).j0 = .5*1089;
            cam(4).i_pitch = .00055;
            cam(4).j_pitch = .00055;
            cam(4).theta = 0;
            cam(4).phi = 0;
            cam(4).bf = 10;
            cam(4).lens_type = 1;
            cam(4).win_d = 15*2.54;
            cam(4).win_t = 2*2.54;
            cam(4).win_phi = pi/6;
            cam(4).win_pitch = 0;
            cam(4).win_yaw = pi/6;
            
            jar_phi = 0.0654;%5.5*pi/180;
            liquidlevel = 51.9310;%53;

            jar_wall = .5;
            jar_OD = 30;
            jar_axwall = .5;
            jar_Oaxrad = 30;
            jar_knucklewall = .5;
            jar_knuckleOrad = 5;
            jar_xpitch = 0;
            jar_ypitch = 0;

            n_air = 1;
            n_window = 1.45846;
            n_hydraulic = 1.430;
            n_jar = 1.45846;
            n_target = 1.188528;
            
            z_offset = 18.153;%7.0871;
            
        case 'secondtry'
            cam(1).x = 0.1338;%0;
            cam(1).y = -1.6705;%-2.5;
            cam(1).z = 40.8678;%7.9 + 13*2.54;
            cam(1).pitch = -0.0249;%0;
            cam(1).yaw = -0.5112;%-pi/6;
            cam(1).roll = 1.5981;%pi/2;
            cam(1).f = 0.5526;%.5;
            cam(1).f1 = -0.0255;%-0.02;
            cam(1).f2 = 0;
            cam(1).i0 = .5*2049 - 170;
            cam(1).j0 = .5*1089;
            cam(1).i_pitch = .00055;
            cam(1).j_pitch = .00055;
            cam(1).theta = 0;
            cam(1).phi = 0;
            cam(1).bf = 10;
            cam(1).lens_type = 1;
            cam(1).win_d = 15*2.54;
            cam(1).win_t = 2*2.54;
            cam(1).win_phi = -pi/6;
            cam(1).win_pitch = 0;
            cam(1).win_yaw = -pi/6;

            cam(2).x = 0.0526;%0;
            cam(2).y = -1.5503;%-2.5;
            cam(2).z = 13.0099;%3.5 + 2*2.54;
            cam(2).pitch = 0.0238;%0;
            cam(2).yaw = -0.5266;%-pi/6;
            cam(2).roll = 1.5762;%pi/2;
            cam(2).f = 0.4914;%.5;
            cam(2).f1 = -0.0219;%0;
            cam(2).f2 = 0;
            cam(2).i0 = .5*2049 - 170;
            cam(2).j0 = .5*1089;
            cam(2).i_pitch = .00055;
            cam(2).j_pitch = .00055;
            cam(2).theta = 0;
            cam(2).phi = 0;
            cam(2).bf = 10;
            cam(2).lens_type = 1;
            cam(2).win_d = 15*2.54;
            cam(2).win_t = 2*2.54;
            cam(2).win_phi = -pi/6;
            cam(2).win_pitch = 0;
            cam(2).win_yaw = -pi/6;

            cam(3).x = 0.2420;%0;
            cam(3).y = -2.4496;%-2.5;
            cam(3).z = 40.9291;%3.5 + 13*2.54;
            cam(3).pitch = -0.0058;%0;
            cam(3).yaw = 0.5167;%pi/6;
            cam(3).roll = 1.5652;%pi/2;
            cam(3).f = 0.5389;%.5;
            cam(3).f1 = -0.0217;%0;
            cam(3).f2 = 0;
            cam(3).i0 = .5*2049 - 170;
            cam(3).j0 = .5*1089;
            cam(3).i_pitch = .00055;
            cam(3).j_pitch = .00055;
            cam(3).theta = 0;
            cam(3).phi = 0;
            cam(3).bf = 10;
            cam(3).lens_type = 1;
            cam(3).win_d = 15*2.54;
            cam(3).win_t = 2*2.54;
            cam(3).win_phi = pi/6;
            cam(3).win_pitch = 0;
            cam(3).win_yaw = pi/6;

            cam(4).x = 0.1520;%0;
            cam(4).y = -3.2227;%-2.5;
            cam(4).z = 12.800;%3.5 + 2*2.54;
            cam(4).pitch = 0.0098;%0;
            cam(4).yaw = 0.5325;%pi/6;
            cam(4).roll = 1.6038;%pi/2;
            cam(4).f = 0.5063;%.5;
            cam(4).f1 = -0.0129;%0;
            cam(4).f2 = 0;
            cam(4).i0 = .5*2049 - 170;
            cam(4).j0 = .5*1089;
            cam(4).i_pitch = .00055;
            cam(4).j_pitch = .00055;
            cam(4).theta = 0;
            cam(4).phi = 0;
            cam(4).bf = 10;
            cam(4).lens_type = 1;
            cam(4).win_d = 15*2.54;
            cam(4).win_t = 2*2.54;
            cam(4).win_phi = pi/6;
            cam(4).win_pitch = 0;
            cam(4).win_yaw = pi/6;
            
            jar_phi = 0.0401;%5.5*pi/180;
            liquidlevel = 52.7941;%53;

            jar_wall = .5;
            jar_OD = 30;
            jar_axwall = .5;
            jar_Oaxrad = 30;
            jar_knucklewall = .5;
            jar_knuckleOrad = 5;
            jar_xpitch = 0;
            jar_ypitch = 0;

            n_air = 1;
            n_window = 1.45846;
            n_hydraulic = 1.430;
            n_jar = 1.45846;
            n_target = 1.188528;
            
            z_offset = 18.0329;%7.0871;
    end

    param_list = [];
    fn = fieldnames(cam);
    for nc=1:length(cam)
        for nf=1:length(fn)
            param_list = [param_list, cam(nc).(fn{nf})];
        end
    end

    paramstruct.(paramsets{np}) = ...
        [param_list, jar_phi, liquidlevel, jar_wall, jar_OD, ...
        jar_axwall, jar_Oaxrad, jar_knucklewall, jar_knuckleOrad, jar_xpitch, jar_ypitch, ...
        n_air, n_window, n_hydraulic, n_jar, n_target, z_offset];

end

if freshfit || ~isfield(paramstruct,'bestfit')
    paramstruct.thisfit = paramstruct.(paramsets{3});
else
    paramstruct.thisfit = paramstruct.bestfit;
end

%% wait, lets add some two-pixel points
% '20160915_1'  5
%%% '20160915_1'        8
% '20160915_1'  11
% '20160915_1'  14
% '20160915_1'  18
%%% '20160915_1'        20
%%% '20160915_1'        23
%%% '20160915_1'        27
%%% '20160915_1'        40
% '20160915_1'  53
% % '20160912_0'        15
% % '20160912_4'        66
%%% % '20160914_1'      16
% % '20160914_1'        40
%%% % '20160912_0'      9
% % '20160912_4'        4
% % '20160912_4'        20
%%% % '20160912_4'      49
%%% % '20160914_1'      28
%%% % '20160914_1'      56
%%% % '20160915_1'      64

% this will be a cell array, each cell an nx4 matrix of images of the same
% wall event.  colums are cam num, ipix, jpix, surface index (2 front, 3
% back)

two_pixel_list = reshape([ ... 20160915_1, 5
        [1, 1060, 597, 3 ; ...
        2, 501, 608, 3 ; ...
        3, 1084, 346, 3 ; ...
        4, 454, 365, 3 ], ...
        ...
        [1, 1395, 755, 2 ; ... 20160915_1, 11
        2, 113, 754, 2 ; ...
        3, 1201, 934, 3 ; ...
        4, 375, 885, 3 ], ...
        ...
        [1, 1047, 309, 3 ; ... 20160915_1, 14
        2, 454, 362, 3 ; ...
        3, 1145, 177, 3 ; ...
        4, 290, 234, 3 ], ...
        ...
        [1, 1270, 643, 3 ; ... 20160915_1, 18
        2, 678, 659, 3 ; ...
        3, 1274, 401, 3 ; ...
        4, 652, 397, 3 ], ...
        ...
        [1, nan, nan, 2 ; ... 20160915_1, 53
        2, 489, 898, 2 ; ...
        3, 1378, 888, 3 ; ...
        4, 627, 866, 3 ], ...
        ...
        [1, 1046, 132, 3 ; ... 20160912_0, 15
         2, 258, 229, 3; ...
         3, 1183, 351, 2 ; ...
         4, nan, nan, 2 ], ...
         ...
         [1, 815, 237, 3 ; ... 20160912_4, 66
         2, 238, 318, 3 ; ...
         3, 843, 148, 3 ; ...
         4, nan, nan, 3 ], ...
         ...
         [1, 1526, 214, 3 ; ... 20160914_1, 40
         2, 821, 237, 3 ; ...
         3, nan, nan, 2 ; ...
         4, 758, 170, 2 ], ...
         ...
         [1, 1420, 693, 2 ; ... 20160912_4, 4
         2, 117, 698, 2 ; ...
         3, 1220, 940, 3 ; ...
         4, 370, 892, 3 ], ...
         ...
         [1, 1073, 962, 3 ; ... 20160912_4, 20
         2, nan, nan, 3 ; ...
         3, 1004, 859, 3 ; ...
         4, 347, 813, 3 ], ...
     ...
     ], 4, 4, []);


twopix_camorder = [1,2,1,3,1,4,2,3,2,4,3,4];

twopix_pix = reshape(permute( ...
    two_pixel_list(twopix_camorder, 1:3, :), ...
    [1, 3, 2]), [], 3);

twopix_izphi = zeros(size(twopix_pix));
twopix_izphi(:, 1) = -1 * reshape( ...
    two_pixel_list(twopix_camorder, 4, :), [], 1);
twopix_izphi(:, 2) = 1;
twopix_izphi(2:2:end, 2) = -1;
twopix_izphi(:, 3) = inf;

twopix_err = zeros(size(twopix_izphi, 1), 1) + 0.1;


badpix = any(isnan(twopix_pix), 2);
badpix(1:2:end) = badpix(1:2:end) | badpix(2:2:end);
badpix(2:2:end) = badpix(1:2:end);

twopix_pix = twopix_pix(~badpix, :);
twopix_izphi = twopix_izphi(~badpix, :);
twopix_err = twopix_err(~badpix);

%% ready to rock!
all_stuffs = [ ...
    [fid_izphi, fid_err, fid_pix] ; ...
    ...[apex_izphi, apex_err, apex_pix] ; ...
    [frontsurface_izphi , frontsurface_err , frontsurface_pix] ; ...
    [backsurface_izphi , backsurface_err , backsurface_pix] ; ...
    [odtan_izphi , odtan_err , odtan_pix] ; ...
    ...[idtan_izphi , idtan_err , idtan_pix] ; ...
    ...[crit_izphi , crit_err , crit_pix] ; ...
    [twopix_izphi, twopix_err, twopix_pix] ; ...
    ];

cam_cut = false(size(all_stuffs,1), length(cam));
for i_c=1:length(cam)
    cam_cut(:, i_c) = all_stuffs(:,5)==i_c;
end

whichcams = [1:4];%[1, 2, 3, 4];
all_camcut = any(cam_cut(:, whichcams), 2);

all_izphi = all_stuffs(all_camcut, 1:3);
all_errs = all_stuffs(all_camcut, 4);
all_pixel_list = all_stuffs(all_camcut, 5:7);

whichparams = false(size(paramstruct.thisfit));


fitparams = paramstruct.thisfit(whichparams);

tic;
res = EricsAutomaticGeometryGoodnessWithTorus_lsq(fitparams, paramstruct.thisfit, whichparams, ...
    all_pixel_list, all_izphi, all_errs, 'float');
chisq = sum(res.^2);
toc;

fprintf(1,'Starting at chisq=%f\n',chisq);
% return

%% really try fitting now....
cam_basicparams = [1:9 ; 23:31 ; 45:53 ; 67:75];
jar_basicparams = [89, 90, 104];
refraction_params = 99:103;
% whichparams([1:7,23:29,45:46]) = true;
% whichparams([8, 30, 50]) = true;
% whichparams([cam_basicparams(:)', jar_basicparams([1, 3])]) = true;
% whichparams([cam_basicparams(4, [1:8])]) = true;
whichparams([reshape(cam_basicparams([1:4], [1:8]), 1, []), jar_basicparams([1:3])]) = true;
fitparams = paramstruct.thisfit(whichparams);

params_ub = paramstruct.upperbound(whichparams);
params_lb = paramstruct.lowerbound(whichparams);

% fitfun = @(params)EricsAutomaticGeometryGoodness(params, paramstruct.thisfit, whichparams, ...
%     all_pixel_list, all_izphi, all_errs, 'arc');
fitfun_lsq = @(params)EricsAutomaticGeometryGoodnessWithTorus_lsq(params, paramstruct.thisfit, whichparams, ...
    all_pixel_list, all_izphi, all_errs, 'float');
fitfun = @(params)(sum(fitfun_lsq(params).^2));

chisq = fitfun(fitparams);

options = optimset('MaxFunEvals',1e3*length(fitparams),'MaxIter',500*length(fitparams),'GradObj','off','Hessian','off');

if 0
    bestfit = fminunc(fitfun,fitparams,options);
    paramstruct.bestfit = paramstruct.thisfit;
    paramstruct.bestfit(whichparams) = bestfit;
elseif 1
    bestfit = lsqnonlin(fitfun_lsq,fitparams,params_lb,params_ub,options);
    paramstruct.bestfit = paramstruct.thisfit;
    paramstruct.bestfit(whichparams) = bestfit;
end

disp(fitfun(fitparams)-fitfun(bestfit));



%% make the lookup table!
if 1
    lookupdirlist = { ...
        '/Users/cdahl/Desktop', ...
        '/nashome/c/cdahl', ...
        '/home/cdahl', ...
        };
    
    for i_d=1:length(lookupdirlist)
        if exist(lookupdirlist{i_d}, 'dir')
            break
        end
    end
    
    lookupdir = lookupdirlist{i_d};
    imagedir = [lookupdir filesep '1' filesep 'Images'];
    
    EricsQuickLookupTableMakerEtcWithTorus(paramstruct.bestfit, ...
        repmat([1700,1088],4,1), ...
        [lookupdir filesep 'PICO60SNOLABlookuptable_Oct2016v2.mat'], ...
        {[imagedir filesep 'cam0_image30.png'], ...
        [imagedir filesep 'cam1_image30.png'], ...
        [imagedir filesep 'cam2_image30.png'], ...
        [imagedir filesep 'cam3_image30.png']});
end
 