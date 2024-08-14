float SFBDT_weight_Zmumu(float bdt, int only_year=0){

  if (only_year==16){
    if (bdt<0.3) return 0.982150556647;
    else if (bdt<0.7) return 0.97966070511;
    else if (bdt<0.9) return 0.9976113738;
    else if (bdt<1.01) return 1.03180754007;
    else return 1; } 
  else if (only_year==17){
    if (bdt<0.3) return 1.20923988846;
    else if (bdt<0.7) return 1.02295653678;
    else if (bdt<0.9) return 0.980403928299;
    else if (bdt<1.01) return 0.973208854917;
    else return 1; } 
  else if (only_year==18){
    if (bdt<0.3) return 1.16067266345;
    else if (bdt<0.7) return 1.00453936862;
    else if (bdt<0.9) return 0.980437576124;
    else if (bdt<1.01) return 1.00998235168;
    else return 1; }
  else{
    if (bdt<0.3) return (0.98*36.3+1.21*41.5+1.16*59.8)/137.6;
    else if (bdt<0.7) return (0.98*36.3+1.02*41.5+1.00*59.8)/137.6;
    else if (bdt<0.9) return (1.00*36.3+0.98*41.5+0.98*59.8)/137.6;
    else if (bdt<1.01) return (1.03*36.3+0.97*41.5+1.01*59.8)/137.6;
    else return 1; }  
}



float SFBDT_weight_Zee(float bdt,int only_year=0){

  if (only_year==16){
    if (bdt<0.3) return 1.0289637271;
    else if (bdt<0.7) return 0.990684014053;
    else if (bdt<0.9) return 0.996662865505;
    else if (bdt<1.01) return 1.01207698081;
    else return 1;} 
  else if (only_year==17){
    if (bdt<0.3) return 1.19244458132;
    else if (bdt<0.7) return 1.01767194623;
    else if (bdt<0.9) return 0.981117820616;
    else if (bdt<1.01) return 0.969830766635;
    else return 1;} 
  else if (only_year==18){
    if (bdt<0.3) return 1.14804193773;
    else if (bdt<0.7) return 1.01099644423;
    else if (bdt<0.9) return 0.982824260351;
    else if (bdt<1.01) return 0.996038455209;
    else return 1;}
  else{
    if (bdt<0.3) return (1.02*36.3+1.19*41.5+1.15*59.8)/137.6;
    else if (bdt<0.7) return (0.99*36.3+1.02*41.5+1.01*59.8)/137.6;
    else if (bdt<0.9) return (1.00*36.3+0.98*41.5+0.98*59.8)/137.6;
    else if (bdt<1.01) return (1.01*36.3+0.97*41.5+1.00*59.8)/137.6;
    else return 1; }  
 }
