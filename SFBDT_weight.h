float SFBDT_weight_Zmumu(float bdt){
if (bdt<0.3) return 1.01693505124;
else if (bdt<0.7) return 0.920792402325;
else if (bdt<0.9) return 1.00877960594;
else if (bdt<1.01) return 1.06255554041;
else return 1;
}
float SFBDT_weight_Zee(float bdt){
if (bdt<0.3) return 1.00646616173;
else if (bdt<0.7) return 0.92343599443;
else if (bdt<0.9) return 1.01730971669;
else if (bdt<1.01) return 1.06366933971;
else return 1;
}
float SFBDT_weightUp_Zmumu(float bdt){
if (bdt<0.3) return 1.02132493854;
else if (bdt<0.7) return 0.922395279172;
else if (bdt<0.9) return 1.00977849261;
else if (bdt<1.01) return 1.06456667338;
else return 1;
}
float SFBDT_weightUp_Zee(float bdt){
if (bdt<0.3) return 1.01079244504;
else if (bdt<0.7) return 0.925087607668;
else if (bdt<0.9) return 1.0183720833;
else if (bdt<1.01) return 1.06597850229;
else return 1;
}
float SFBDT_weightDown_Zmumu(float bdt){
if (bdt<0.3) return 1.01254516393;
else if (bdt<0.7) return 0.919189525478;
else if (bdt<0.9) return 1.00778071928;
else if (bdt<1.01) return 1.06054440743;
else return 1;
}
float SFBDT_weightDown_Zee(float bdt){
if (bdt<0.3) return 1.00213987842;
else if (bdt<0.7) return 0.921784381192;
else if (bdt<0.9) return 1.01624735009;
else if (bdt<1.01) return 1.06136017714;
else return 1;
}
