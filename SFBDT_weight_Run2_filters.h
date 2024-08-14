float SFBDT_weight_Zmumu(float bdt){
if (bdt<0.3) return 1.11976629482;
else if (bdt<0.7) return 1.0007865426;
else if (bdt<0.9) return 0.985680889785;
else if (bdt<1.01) return 1.00754063889;
else return 1;
}
float SFBDT_weight_Zee(float bdt){
if (bdt<0.3) return 1.11685522174;
else if (bdt<0.7) return 1.00175273685;
else if (bdt<0.9) return 0.985780358772;
else if (bdt<1.01) return 0.998934559096;
else return 1;
}

float SFBDT_weightUp_Zmumu(float bdt){
if (bdt<0.3) return 1.12334312979;
else if (bdt<0.7) return 1.00207600662;
else if (bdt<0.9) return 0.986373148003;
else if (bdt<1.01) return 1.00890201363;
else return 1;
}
float SFBDT_weightUp_Zee(float bdt){
if (bdt<0.3) return 1.12057700491;
else if (bdt<0.7) return 1.00313273497;
else if (bdt<0.9) return 0.98654583458;
else if (bdt<1.01) return 1.00054308706;
else return 1;
}
float SFBDT_weightDown_Zmumu(float bdt){
if (bdt<0.3) return 1.11618945985;
else if (bdt<0.7) return 0.99949707858;
else if (bdt<0.9) return 0.984988631567;
else if (bdt<1.01) return 1.00617926414;
else return 1;
}
float SFBDT_weightDown_Zee(float bdt){
if (bdt<0.3) return 1.11313343858;
else if (bdt<0.7) return 1.00037273872;
else if (bdt<0.9) return 0.985014882965;
else if (bdt<1.01) return 0.997326031131;
else return 1;
}
