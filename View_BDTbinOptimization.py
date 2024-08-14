import ROOT as rt

result = "result_combine_4bin_0p2.txt"

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


bdt_cuts=[]
limits=[]
with open(result,"r") as txt:
  lines = txt.readlines()
  temp_bdt=[]
  for line in lines:
    temp_bdt=[]
    words = line.split(" ")
    if "bdt_cut" in words:
       for word in words:
          if is_number(word):
             temp_bdt.append( float(word) )
       bdt_cuts.append(temp_bdt)
    if "Observed" in words:
       limits.append(float( words[-1]))
   
for idx in range(len(bdt_cuts)):
  print bdt_cuts[idx],limits[idx]
    

print "min",limits[limits.index(min(limits))],"cuts",bdt_cuts[limits.index(min(limits))]

