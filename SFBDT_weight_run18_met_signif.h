float SFBDT_weight_Zmumu(float bdt){
if (bdt<0.3) return 1.36465487879;
else if (bdt<0.7) return 1.0348608341;
else if (bdt<0.9) return 0.965826255541;
else if (bdt<1.01) return 0.985168164657;
else return 1;
}
float SFBDT_weight_Zee(float bdt){
if (bdt<0.3) return 1.37266334693;
else if (bdt<0.7) return 1.03750841915;
else if (bdt<0.9) return 0.966097047613;
else if (bdt<1.01) return 0.968319390266;
else return 1;
}
float SFBDT_weight_Zll(float bdt){
if (bdt<0.3) return (1.37266334693+1.36465487879)/2.;
else if (bdt<0.7) return (1.03750841915+1.0348608341)/2.;
else if (bdt<0.9) return (0.966097047613+0.965826255541)/2.;
else if (bdt<1.01) return (0.968319390266+0.985168164657)/2.;
else return 1;
}
float SFBDT_weightUp_Zmumu(float bdt){
if (bdt<0.3) return 1.37219613446;
else if (bdt<0.7) return 1.03701502583;
else if (bdt<0.9) return 0.966875103642;
else if (bdt<1.01) return 0.987209728801;
else return 1;
}
float SFBDT_weightUp_Zee(float bdt){
if (bdt<0.3) return 1.38011679066;
else if (bdt<0.7) return 1.03970369929;
else if (bdt<0.9) return 0.967208584187;
else if (bdt<1.01) return 0.970631993901;
else return 1;
}
float SFBDT_weightDown_Zmumu(float bdt){
if (bdt<0.3) return 1.35711362312;
else if (bdt<0.7) return 1.03270664237;
else if (bdt<0.9) return 0.96477740744;
else if (bdt<1.01) return 0.983126600513;
else return 1;
}
float SFBDT_weightDown_Zee(float bdt){
if (bdt<0.3) return 1.3652099032;
else if (bdt<0.7) return 1.035313139;
else if (bdt<0.9) return 0.964985511038;
else if (bdt<1.01) return 0.966006786631;
else return 1;
}
