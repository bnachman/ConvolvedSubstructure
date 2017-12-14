from ROOT import *
import numpy as np

#We are going to move a high mass D2 distribution to a low mass D2 distribution.  By construction, the low mass distribution should not change.

myfile = TFile("myfile.root")
mytree = myfile.Get("mytree")
D2histo_lowm = TH1D("low_m","low_m",50,0,5)
D2histo_highm = TH1D("high_m","high_m",50,0,5)

lowm = 50.
highm = 200.
epsilon = 10.

for i in range(mytree.GetEntries()):
  if (mytree.m > lowm-epsilon and mytree.m < lowm+epsilon):
    D2histo_lowm.Fill(mytree.D2)
  elif (mytree.m > highm-epsilon and mytree.m < highm+epsilon):
    D2histo_highm.Fill(mytree.D2)
    pass
  pass

#This is the DDT version
D2histo_lowm_ddt = TH1D("low_m_ddt","low_m_ddt",50,0,5)
D2histo_highm_ddt = TH1D("high_m_ddt","high_m_ddt",50,0,5)
for i in range(mytree.GetEntries()):
  if (mytree.m > lowm-epsilon and mytree.m < lowm+epsilon):
    D2histo_lowm_ddt.Fill(mytree.D2-D2histo_lowm.GetMean()+D2histo_lowm.GetMean())
  elif (mytree.m > highm-epsilon and mytree.m < highm+epsilon):
    D2histo_highm_ddt.Fill(mytree.D2-D2histo_highm.GetMean()+D2histo_lowm.GetMean())
    pass
  pass

shapefunc_low = D2histo_lowm.Clone()
shapefunc_high = D2histo_highm.Clone()

#These are free parameters - need to tune these.
Omega_low = 0.01 
Omega_high = 0.23
shapeval = 2.4

#First, let's make the shape functions that we are going to use to convolve with the D2 distributions.
for i in range(1,shapefunc_low.GetNbinsX()+1):
  xx = shapefunc_low.GetXaxis().GetBinCenter(i)

  myval = (shapeval**shapeval/ROOT.Math.tgamma(shapeval))*(xx**(shapeval-1))*np.exp(-shapeval*xx/Omega_low) / (Omega_low)**shapeval
  shapefunc_low.SetBinContent(i,myval)

  myval = (shapeval**shapeval/ROOT.Math.tgamma(shapeval))*(xx**(shapeval-1))*np.exp(-shapeval*xx/Omega_high) / (Omega_high)**shapeval
  shapefunc_high.SetBinContent(i,myval)

  pass

print "Sanity check! These should be 1: ",shapefunc_low.Integral("width"),shapefunc_high.Integral("width")

#Now, let's do some convolving!
D2histo_lowm_conv = D2histo_lowm.Clone("D2histo_lowm_conv") 
D2histo_highm_conv = D2histo_highm.Clone("D2histo_highm_conv")  
for i in range(0,shapefunc_low.GetNbinsX()):
  mysum = 0.
  mysum2 = 0.
  for j in range(0,shapefunc_low.GetNbinsX()):
    if (i-j >= 0):
      mysum+=D2histo_lowm_conv.GetBinContent(j+1)*shapefunc_low.GetBinContent(i-j+1)
      mysum2+=D2histo_highm_conv.GetBinContent(j+1)*shapefunc_high.GetBinContent(i-j+1)
      pass
    pass
  D2histo_lowm_conv.SetBinContent(i+1,D2histo_lowm.GetBinContent(i+1))
  #D2histo_lowm_conv.SetBinContent(i,mysum)
  D2histo_highm_conv.SetBinContent(i+1,mysum2)
  pass

#Now, let's the actual transformation for D2, i.e. D2 -> D2_CSS
D2histo_lowm_CDF = D2histo_lowm.Clone("D2histo_lowm_CDF")
D2histo_lowm_conv_CDF = D2histo_lowm_conv.Clone("D2histo_lowm_conv")

D2histo_highm_CDF = D2histo_highm.Clone("D2histo_highm_CDF")
D2histo_highm_conv_CDF = D2histo_highm_conv.Clone("D2histo_highm_conv")

for i in range(1,D2histo_lowm_CDF.GetNbinsX()+1):
  D2histo_lowm_CDF.SetBinContent(i,D2histo_lowm.Integral(1,i))
  D2histo_lowm_conv_CDF.SetBinContent(i,D2histo_lowm_conv.Integral(1,i))

  D2histo_highm_CDF.SetBinContent(i,D2histo_highm.Integral(1,i))
  D2histo_highm_conv_CDF.SetBinContent(i,D2histo_highm_conv.Integral(1,i))   
  pass

#Let F = CDF of original function and G = CDF of convolution
#Then, our map shoud be G^{-1}(F(X)).
#http://math.arizona.edu/~jwatkins/f-transform.pdf page 3.

fxvals = []
fyvals = []
ginvxvals = []
ginvyvals = []

fxvals2 = []
fyvals2 = []
ginvxvals2 = []
ginvyvals2 = []

for i in range(1,D2histo_lowm_CDF.GetNbinsX()+1):
  xx = D2histo_lowm_CDF.GetXaxis().GetBinCenter(i)
  fxvals+=[xx]
  ginvyvals+=[xx]
  fyvals+=[D2histo_lowm_CDF.GetBinContent(i)]
  ginvxvals+=[D2histo_lowm_conv_CDF.GetBinContent(i)]

  fxvals2+=[xx]
  ginvyvals2+=[xx]
  fyvals2+=[D2histo_highm_CDF.GetBinContent(i)]
  ginvxvals2+=[D2histo_highm_conv_CDF.GetBinContent(i)]

  pass

F = TGraph(len(fxvals),np.array(fxvals),np.array(fyvals))
Ginv = TGraph(len(ginvxvals),np.array(ginvxvals),np.array(ginvyvals))

F2 = TGraph(len(fxvals2),np.array(fxvals2),np.array(fyvals2))
Ginv2 = TGraph(len(ginvxvals2),np.array(ginvxvals2),np.array(ginvyvals2))

#This is the CSS version
D2histo_lowm_css = TH1D("low_m_css","low_m_css",50,0,5)
D2histo_highm_css = TH1D("high_m_css","high_m_css",50,0,5)
for i in range(mytree.GetEntries()):
  if (mytree.m > lowm-epsilon and mytree.m < lowm+epsilon):
    D2histo_lowm_css.Fill(Ginv.Eval(F.Eval(mytree.D2)))
  elif (mytree.m > highm-epsilon and mytree.m < highm+epsilon):
    D2histo_highm_css.Fill(Ginv2.Eval(F2.Eval(mytree.D2)))
    pass
  pass
