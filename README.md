# nlmp_bvp

**Nonlinear Multipoint Boundary Value Problem Solver** in `C++` based on:  
1. Ojika, T., and Y. Kasue. "Initial-value adjusting method for the solution of nonlinear multipoint boundary-value problems."  
Journal of Mathematical Analysis and Applications 69.2 (1979): 359-371.  
2. Welsh, Wayne, and Takeo Ojika. "Multipoint boundary value problems with discontinuities I. Algorithms and applications."  
Journal of Computational and Applied Mathematics 6.2 (1980): 133-143.  
This was in response to a need for a suitable C++ alternative to bvp4c in Mathworks MATLAB.  

## Output - test_nlmpBVP.cpp:
```
[shivanandvp@computer nlmpbvp]$ /mnt/2E42F6A042F66BC9/Development/nlmpbvp/build/test_nlmpBVP

=============================================================================
Test: Non-linear multipoint boundary value problem solver (nlmpBVP, nlmpBVP2)
=============================================================================
Copyright shivanandvp (shivanandvp.oss@gmail.com)

========================
Boundary Value Problem 1
========================

Boundary nodes (tBC) = 
0 4

Starting state vector (oxt1) = 
1
0

Initiating the BVP solver...

Boundary nodes correspond to the below columns: 
  0 100

k = 0... kG = 2.66286... kGPrev = inf alpha = 1 Increasing alpha...
k = 1... kG = 5.24021... kGPrev = 2.66286 alpha = 1 Decreasing alpha...
k = 2... kG = 1.81289e-08... kGPrev = 5.24021 alpha = 1 Increasing alpha...
Ran 3 iteration(s).

========
Solution
========
Independent variable (t) = 
***REMOVED***0 0.04 0.08 0.12 0.16  0.2 0.24 0.28 0.32 0.36  0.4 0.44 0.48 0.52 0.56  0.6 0.64 0.68 0.72 0.76  0.8 0.84 0.88 0.92 0.96***REMOVED*** 1 1.04 1.08 1.12 1.16  1.2 1.24 1.28 1.32 1.36  1.4 1.44 1.48 1.52 1.56  1.6 1.64 1.68 1.72 1.76  1.8 1.84 1.88 1.92 1.96***REMOVED*** 2 2.04 2.08 2.12 2.16  2.2 2.24 2.28 2.32 2.36  2.4 2.44 2.48 2.52 2.56  2.6 2.64 2.68 2.72 2.76  2.8 2.84 2.88 2.92 2.96***REMOVED*** 3 3.04 3.08 3.12 3.16  3.2 3.24 3.28 3.32 3.36  3.4 3.44 3.48 3.52 3.56  3.6 3.64 3.68 3.72 3.76  3.8 3.84 3.88 3.92 3.96***REMOVED*** 4

State vector (x) = 
***REMOVED******REMOVED******REMOVED*** 0 -0.00293227 -0.00586923 -0.00881558***REMOVED***-0.011776  -0.0147553  -0.0177583  -0.0207896  -0.0238542***REMOVED***-0.026957  -0.0301029  -0.0332969  -0.0365443  -0.0398501  -0.0432198  -0.0466585***REMOVED***-0.050172  -0.0537657  -0.0574454  -0.0612171  -0.0650867  -0.0690605  -0.0731448  -0.0773462  -0.0816713  -0.0861271  -0.0907208  -0.0954596***REMOVED***-0.100351***REMOVED***-0.105403***REMOVED***-0.110624***REMOVED***-0.116022***REMOVED***-0.121605***REMOVED***-0.127384***REMOVED***-0.133365***REMOVED***-0.139561***REMOVED***-0.145979***REMOVED***-0.152632***REMOVED***-0.159528***REMOVED*** -0.16668***REMOVED***-0.174099***REMOVED***-0.181796***REMOVED***-0.189784***REMOVED***-0.198075***REMOVED***-0.206684***REMOVED***-0.215624***REMOVED***-0.224908***REMOVED***-0.234552***REMOVED***-0.244572***REMOVED***-0.254983***REMOVED***-0.265802***REMOVED***-0.277047***REMOVED***-0.288734***REMOVED***-0.300884***REMOVED***-0.313516***REMOVED***-0.326648***REMOVED***-0.340304***REMOVED***-0.354504***REMOVED***-0.369272***REMOVED*** -0.38463***REMOVED***-0.400604***REMOVED***-0.417219***REMOVED***-0.434502***REMOVED*** -0.45248***REMOVED***-0.471182***REMOVED***-0.490638***REMOVED***-0.510879***REMOVED***-0.531937***REMOVED***-0.553847***REMOVED***-0.576643***REMOVED***-0.600362***REMOVED***-0.625042***REMOVED***-0.650722***REMOVED***-0.677443***REMOVED***-0.705248***REMOVED***-0.734181***REMOVED*** -0.76429***REMOVED***-0.795621***REMOVED***-0.828226***REMOVED***-0.862156***REMOVED***-0.897466***REMOVED***-0.934212***REMOVED***-0.972452***REMOVED*** -1.01225***REMOVED*** -1.05367***REMOVED*** -1.09677***REMOVED*** -1.14163***REMOVED*** -1.18831***REMOVED***  -1.2369***REMOVED*** -1.28746***REMOVED*** -1.34009***REMOVED*** -1.39486***REMOVED*** -1.45186***REMOVED*** -1.51119***REMOVED*** -1.57293***REMOVED*** -1.63719***REMOVED*** -1.70407***REMOVED*** -1.77368***REMOVED*** -1.84613***REMOVED*** -1.92153***REMOVED******REMOVED******REMOVED*** -2
 -0.0732871  -0.0733458  -0.0735218  -0.0738154  -0.0742272  -0.0747578***REMOVED***-0.075408  -0.0761788  -0.0770716  -0.0780877  -0.0792287  -0.0804965  -0.0818932  -0.0834209***REMOVED***-0.085082  -0.0868794  -0.0888157  -0.0908942  -0.0931181***REMOVED***-0.095491  -0.0980168***REMOVED***-0.100699***REMOVED***-0.103543***REMOVED***-0.106553***REMOVED***-0.109732***REMOVED***-0.113088***REMOVED***-0.116624***REMOVED***-0.120348***REMOVED***-0.124263***REMOVED***-0.128378***REMOVED***-0.132698***REMOVED*** -0.13723***REMOVED***-0.141982***REMOVED***-0.146961***REMOVED***-0.152175***REMOVED***-0.157633***REMOVED***-0.163343***REMOVED***-0.169315***REMOVED***-0.175557***REMOVED*** -0.18208***REMOVED***-0.188895***REMOVED***-0.196012***REMOVED***-0.203443***REMOVED***-0.211199***REMOVED***-0.219293***REMOVED***-0.227738***REMOVED***-0.236547***REMOVED***-0.245735***REMOVED***-0.255316***REMOVED***-0.265306***REMOVED***-0.275721***REMOVED***-0.286576***REMOVED*** -0.29789***REMOVED***-0.309681***REMOVED***-0.321967***REMOVED***-0.334769***REMOVED***-0.348106***REMOVED***-0.362001***REMOVED***-0.376474***REMOVED*** -0.39155***REMOVED***-0.407253***REMOVED***-0.423607***REMOVED***-0.440639***REMOVED***-0.458376***REMOVED***-0.476847***REMOVED***-0.496081***REMOVED***-0.516109***REMOVED***-0.536962***REMOVED***-0.558675***REMOVED***-0.581282***REMOVED***-0.604819***REMOVED***-0.629324***REMOVED***-0.654836***REMOVED***-0.681395***REMOVED***-0.709045***REMOVED*** -0.73783***REMOVED***-0.767796***REMOVED*** -0.79899***REMOVED***-0.831462***REMOVED***-0.865265***REMOVED***-0.900453***REMOVED***-0.937082***REMOVED*** -0.97521***REMOVED***  -1.0149***REMOVED*** -1.05621***REMOVED*** -1.09921***REMOVED*** -1.14398***REMOVED*** -1.19057***REMOVED*** -1.23907***REMOVED*** -1.28955***REMOVED*** -1.34209***REMOVED*** -1.39678***REMOVED*** -1.45371***REMOVED*** -1.51296***REMOVED*** -1.57464***REMOVED*** -1.63883***REMOVED*** -1.70565***REMOVED*** -1.77519***REMOVED*** -1.84758***REMOVED*** -1.92292***REMOVED*** -2.00134

Independent variable at boundary nodes (tBC) = 0 4

State vector at boundary nodes (xBC) = 
***REMOVED******REMOVED******REMOVED***0***REMOVED******REMOVED******REMOVED***-2
-0.0732871***REMOVED***-2.00134

====================================================================================================

========================
Boundary Value Problem 2
========================

Boundary nodes (tBC) = 
  0 0.5***REMOVED***1

Starting state vector (oxt1) = 
1
1
1

Initiating the BVP solver...

Boundary nodes correspond to the below columns: 
  0  50 100

k = 0... kG = 49.7451... kGPrev = inf alpha = 1 Increasing alpha...
k = 1... kG = 3.51138e-07... kGPrev = 49.7451 alpha = 1 Increasing alpha...
Ran 2 iteration(s).

========
Solution
========
Independent variable (t) = 
***REMOVED***0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09  0.1 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19  0.2 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29  0.3 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39  0.4 0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49  0.5 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59  0.6 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69  0.7 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79  0.8 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89  0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99***REMOVED*** 1

State vector (x) = 
  -0.0121071***REMOVED***-0.0120974***REMOVED***-0.0120689***REMOVED***-0.0120226***REMOVED***-0.0119594***REMOVED***  -0.01188***REMOVED***-0.0117854***REMOVED***-0.0116762***REMOVED***-0.0115532***REMOVED***-0.0114171***REMOVED***-0.0112685***REMOVED***-0.0111081***REMOVED***-0.0109365***REMOVED***-0.0107541***REMOVED***-0.0105617***REMOVED***-0.0103597***REMOVED***-0.0101486  -0.00992884  -0.00970091  -0.00946523  -0.00922221  -0.00897224  -0.00871571  -0.00845296  -0.00818434  -0.00791018***REMOVED***-0.0076308***REMOVED***-0.0073465  -0.00705755  -0.00676426  -0.00646687  -0.00616565  -0.00586084  -0.00555268***REMOVED***-0.0052414  -0.00492723  -0.00461037  -0.00429104  -0.00396944  -0.00364576  -0.00332019  -0.00299293  -0.00266414  -0.00233402  -0.00200273  -0.00167045  -0.00133734  -0.00100358 -0.000669325 -0.000334744 -8.77731e-19  0.000334744  0.000669325***REMOVED***0.00100358***REMOVED***0.00133734***REMOVED***0.00167045***REMOVED***0.00200273***REMOVED***0.00233402***REMOVED***0.00266414***REMOVED***0.00299293***REMOVED***0.00332019***REMOVED***0.00364576***REMOVED***0.00396944***REMOVED***0.00429104***REMOVED***0.00461037***REMOVED***0.00492723***REMOVED*** 0.0052414***REMOVED***0.00555268***REMOVED***0.00586084***REMOVED***0.00616565***REMOVED***0.00646687***REMOVED***0.00676426***REMOVED***0.00705755***REMOVED*** 0.0073465***REMOVED*** 0.0076308***REMOVED***0.00791018***REMOVED***0.00818434***REMOVED***0.00845296***REMOVED***0.00871571***REMOVED***0.00897224***REMOVED***0.00922221***REMOVED***0.00946523***REMOVED***0.00970091***REMOVED***0.00992884***REMOVED*** 0.0101486***REMOVED*** 0.0103597***REMOVED*** 0.0105617***REMOVED*** 0.0107541***REMOVED*** 0.0109365***REMOVED*** 0.0111081***REMOVED*** 0.0112685***REMOVED*** 0.0114171***REMOVED*** 0.0115532***REMOVED*** 0.0116762***REMOVED*** 0.0117854***REMOVED******REMOVED***0.01188***REMOVED*** 0.0119594***REMOVED*** 0.0120226***REMOVED*** 0.0120689***REMOVED*** 0.0120974***REMOVED*** 0.0121071
-4.47821e-26***REMOVED***0.00192404***REMOVED***0.00375287***REMOVED***0.00549107***REMOVED***0.00714297***REMOVED***0.00871271***REMOVED*** 0.0102042***REMOVED*** 0.0116212***REMOVED*** 0.0129673***REMOVED*** 0.0142457***REMOVED*** 0.0154598***REMOVED*** 0.0166124***REMOVED*** 0.0177067***REMOVED*** 0.0187451***REMOVED*** 0.0197304***REMOVED***  0.020665***REMOVED*** 0.0215513***REMOVED*** 0.0223915***REMOVED*** 0.0231876***REMOVED*** 0.0239417***REMOVED*** 0.0246556***REMOVED*** 0.0253311***REMOVED******REMOVED***0.02597***REMOVED*** 0.0265738***REMOVED***  0.027144***REMOVED*** 0.0276821***REMOVED*** 0.0281894***REMOVED*** 0.0286671***REMOVED*** 0.0291165***REMOVED*** 0.0295387***REMOVED*** 0.0299347***REMOVED*** 0.0303056***REMOVED*** 0.0306522***REMOVED*** 0.0309754***REMOVED*** 0.0312761***REMOVED***  0.031555***REMOVED*** 0.0318127***REMOVED******REMOVED***0.03205***REMOVED*** 0.0322674***REMOVED*** 0.0324654***REMOVED*** 0.0326447***REMOVED*** 0.0328055***REMOVED*** 0.0329483***REMOVED*** 0.0330735***REMOVED*** 0.0331814***REMOVED*** 0.0332722***REMOVED*** 0.0333463***REMOVED*** 0.0334036***REMOVED*** 0.0334445***REMOVED***  0.033469***REMOVED*** 0.0334772***REMOVED***  0.033469***REMOVED*** 0.0334445***REMOVED*** 0.0334036***REMOVED*** 0.0333463***REMOVED*** 0.0332722***REMOVED*** 0.0331814***REMOVED*** 0.0330735***REMOVED*** 0.0329483***REMOVED*** 0.0328055***REMOVED*** 0.0326447***REMOVED*** 0.0324654***REMOVED*** 0.0322674***REMOVED******REMOVED***0.03205***REMOVED*** 0.0318127***REMOVED***  0.031555***REMOVED*** 0.0312761***REMOVED*** 0.0309754***REMOVED*** 0.0306522***REMOVED*** 0.0303056***REMOVED*** 0.0299347***REMOVED*** 0.0295387***REMOVED*** 0.0291165***REMOVED*** 0.0286671***REMOVED*** 0.0281894***REMOVED*** 0.0276821***REMOVED***  0.027144***REMOVED*** 0.0265738***REMOVED******REMOVED***0.02597***REMOVED*** 0.0253311***REMOVED*** 0.0246556***REMOVED*** 0.0239417***REMOVED*** 0.0231876***REMOVED*** 0.0223915***REMOVED*** 0.0215513***REMOVED***  0.020665***REMOVED*** 0.0197304***REMOVED*** 0.0187451***REMOVED*** 0.0177067***REMOVED*** 0.0166124***REMOVED*** 0.0154598***REMOVED*** 0.0142457***REMOVED*** 0.0129673***REMOVED*** 0.0116212***REMOVED*** 0.0102042***REMOVED***0.00871271***REMOVED***0.00714297***REMOVED***0.00549107***REMOVED***0.00375287***REMOVED***0.00192404  5.78964e-17
***REMOVED*** 0.197323***REMOVED***  0.187565***REMOVED***  0.178277***REMOVED***  0.169434***REMOVED***  0.161015***REMOVED***  0.152999***REMOVED***  0.145365***REMOVED***  0.138095***REMOVED******REMOVED***0.13117***REMOVED***  0.124573***REMOVED***  0.118287***REMOVED***  0.112298***REMOVED***  0.106589***REMOVED***  0.101146***REMOVED*** 0.0959568***REMOVED*** 0.0910073***REMOVED*** 0.0862853***REMOVED*** 0.0817791***REMOVED*** 0.0774774***REMOVED*** 0.0733694***REMOVED*** 0.0694448***REMOVED***  0.065694***REMOVED*** 0.0621074***REMOVED*** 0.0586761***REMOVED*** 0.0553915***REMOVED*** 0.0522454***REMOVED*** 0.0492299***REMOVED*** 0.0463376***REMOVED*** 0.0435611***REMOVED*** 0.0408936***REMOVED*** 0.0383283***REMOVED*** 0.0358589***REMOVED*** 0.0334791***REMOVED***  0.031183***REMOVED*** 0.0289649***REMOVED*** 0.0268192***REMOVED*** 0.0247406***REMOVED*** 0.0227239***REMOVED***  0.020764***REMOVED***  0.018856***REMOVED*** 0.0169951***REMOVED*** 0.0151768***REMOVED*** 0.0133964***REMOVED*** 0.0116495***REMOVED*** 0.0099317***REMOVED***0.00823876***REMOVED***0.00656642***REMOVED*** 0.0049105***REMOVED***0.00326686***REMOVED***0.00163139  7.19398e-12  -0.00163139  -0.00326686***REMOVED***-0.0049105  -0.00656642  -0.00823876***REMOVED***-0.0099317***REMOVED***-0.0116495***REMOVED***-0.0133964***REMOVED***-0.0151768***REMOVED***-0.0169951***REMOVED*** -0.018856***REMOVED*** -0.020764***REMOVED***-0.0227239***REMOVED***-0.0247406***REMOVED***-0.0268192***REMOVED***-0.0289649***REMOVED*** -0.031183***REMOVED***-0.0334791***REMOVED***-0.0358589***REMOVED***-0.0383283***REMOVED***-0.0408936***REMOVED***-0.0435611***REMOVED***-0.0463376***REMOVED***-0.0492299***REMOVED***-0.0522454***REMOVED***-0.0553915***REMOVED***-0.0586761***REMOVED***-0.0621074***REMOVED*** -0.065694***REMOVED***-0.0694448***REMOVED***-0.0733694***REMOVED***-0.0774774***REMOVED***-0.0817791***REMOVED***-0.0862853***REMOVED***-0.0910073***REMOVED***-0.0959568***REMOVED*** -0.101146***REMOVED*** -0.106589***REMOVED*** -0.112298***REMOVED*** -0.118287***REMOVED*** -0.124573***REMOVED***  -0.13117***REMOVED*** -0.138095***REMOVED*** -0.145365***REMOVED*** -0.152999***REMOVED*** -0.161015***REMOVED*** -0.169434***REMOVED*** -0.178277***REMOVED*** -0.187565***REMOVED*** -0.197323

Independent variable at boundary nodes (tBC) =***REMOVED***0 0.5***REMOVED***1

State vector at boundary nodes (xBC) = 
  -0.0121071 -8.77731e-19***REMOVED*** 0.0121071
-4.47821e-26***REMOVED*** 0.0334772  5.78964e-17
***REMOVED*** 0.197323  7.19398e-12***REMOVED*** -0.197323

====================================================================================================

========================
Boundary Value Problem 3
========================

Boundary nodes (tBC) = 
***REMOVED******REMOVED*** 0 0.523599***REMOVED***1.0472***REMOVED***1.5708***REMOVED***2.0944  3.14159

Starting state vector (oxt1) = 
 0.1  0.1  0.4  0.8  0.9
-0.6  0.1  0.9  2.1  0.8

Initiating the BVP solver...

Boundary nodes correspond to the below columns: 
  0  20  40  60  80 120

k = 0... kG = 1.09572... kGPrev = inf alpha = 1 Increasing alpha...
k = 1... kG = 2.92486e-09... kGPrev = 1.09572 alpha = 1 Increasing alpha...
Ran 2 iteration(s).

========
Solution
========
Independent variable (t) = 
***REMOVED******REMOVED***  0 0.0261799 0.0523599 0.0785398***REMOVED***0.10472***REMOVED*** 0.1309***REMOVED***0.15708***REMOVED***0.18326***REMOVED***0.20944  0.235619  0.261799  0.287979  0.314159  0.340339  0.366519  0.392699  0.418879  0.445059  0.471239  0.497419  0.523599  0.523599  0.549779  0.575959  0.602139  0.628319  0.654498  0.680678  0.706858  0.733038  0.759218  0.785398  0.811578  0.837758  0.863938  0.890118  0.916298  0.942478  0.968658  0.994838***REMOVED***1.02102***REMOVED*** 1.0472***REMOVED*** 1.0472***REMOVED***1.07338***REMOVED***1.09956***REMOVED***1.12574***REMOVED***1.15192***REMOVED*** 1.1781***REMOVED***1.20428***REMOVED***1.23046***REMOVED***1.25664***REMOVED***1.28282***REMOVED***  1.309***REMOVED***1.33518***REMOVED***1.36136***REMOVED***1.38754***REMOVED***1.41372***REMOVED*** 1.4399***REMOVED***1.46608***REMOVED***1.49226***REMOVED***1.51844***REMOVED***1.54462***REMOVED*** 1.5708***REMOVED*** 1.5708***REMOVED***1.59698***REMOVED***1.62316***REMOVED***1.64934***REMOVED***1.67552***REMOVED*** 1.7017***REMOVED***1.72788***REMOVED***1.75406***REMOVED***1.78024***REMOVED***1.80642***REMOVED*** 1.8326***REMOVED***1.85878***REMOVED***1.88496***REMOVED***1.91114***REMOVED***1.93732***REMOVED*** 1.9635***REMOVED***1.98968***REMOVED***2.01586***REMOVED***2.04204***REMOVED***2.06822***REMOVED*** 2.0944***REMOVED*** 2.0944***REMOVED***2.12058***REMOVED***2.14675***REMOVED***2.17293***REMOVED***2.19911***REMOVED***2.22529***REMOVED***2.25147***REMOVED***2.27765***REMOVED***2.30383***REMOVED***2.33001***REMOVED***2.35619***REMOVED***2.38237***REMOVED***2.40855***REMOVED***2.43473***REMOVED***2.46091***REMOVED***2.48709***REMOVED***2.51327***REMOVED***2.53945***REMOVED***2.56563***REMOVED***2.59181***REMOVED***2.61799***REMOVED***2.64417***REMOVED***2.67035***REMOVED***2.69653***REMOVED***2.72271***REMOVED***2.74889***REMOVED***2.77507***REMOVED***2.80125***REMOVED***2.82743***REMOVED***2.85361***REMOVED***2.87979***REMOVED***2.90597***REMOVED***2.93215***REMOVED***2.95833***REMOVED***2.98451***REMOVED***3.01069***REMOVED***3.03687***REMOVED***3.06305***REMOVED***3.08923***REMOVED***3.11541***REMOVED***3.14159

State vector (x) = 
 1.40454e-28***REMOVED*** 0.0261769***REMOVED***  0.052336***REMOVED*** 0.0784591***REMOVED***  0.104528***REMOVED***  0.130526***REMOVED***  0.156434***REMOVED***  0.182236***REMOVED***  0.207912***REMOVED***  0.233445***REMOVED***  0.258819***REMOVED***  0.284015***REMOVED***  0.309017***REMOVED***  0.333807***REMOVED***  0.358368***REMOVED***  0.382683***REMOVED***  0.406737***REMOVED***  0.430511***REMOVED******REMOVED***0.45399***REMOVED***  0.477159***REMOVED******REMOVED******REMOVED*** 0.5***REMOVED******REMOVED******REMOVED***-0.5***REMOVED*** -0.477159***REMOVED***  -0.45399***REMOVED*** -0.430511***REMOVED*** -0.406737***REMOVED*** -0.382683***REMOVED*** -0.358368***REMOVED*** -0.333807***REMOVED*** -0.309017***REMOVED*** -0.284015***REMOVED*** -0.258819***REMOVED*** -0.233445***REMOVED*** -0.207912***REMOVED*** -0.182236***REMOVED*** -0.156434***REMOVED*** -0.130526***REMOVED*** -0.104528***REMOVED***-0.0784591***REMOVED*** -0.052336***REMOVED***-0.0261769 -1.42748e-12 -1.42748e-12***REMOVED*** 0.0523539***REMOVED***  0.104672***REMOVED***  0.156918***REMOVED***  0.209057***REMOVED***  0.261052***REMOVED***  0.312869***REMOVED***  0.364471***REMOVED***  0.415823***REMOVED***  0.466891***REMOVED***  0.517638***REMOVED***  0.568031***REMOVED***  0.618034***REMOVED***  0.667614***REMOVED***  0.716736***REMOVED***  0.765367***REMOVED***  0.813473***REMOVED***  0.861022***REMOVED***  0.907981***REMOVED***  0.954318***REMOVED******REMOVED******REMOVED******REMOVED***1 -4.25982e-12***REMOVED*** 0.0261769***REMOVED***  0.052336***REMOVED*** 0.0784591***REMOVED***  0.104528***REMOVED***  0.130526***REMOVED***  0.156434***REMOVED***  0.182236***REMOVED***  0.207912***REMOVED***  0.233445***REMOVED***  0.258819***REMOVED***  0.284015***REMOVED***  0.309017***REMOVED***  0.333807***REMOVED***  0.358368***REMOVED***  0.382683***REMOVED***  0.406737***REMOVED***  0.430511***REMOVED******REMOVED***0.45399***REMOVED***  0.477159***REMOVED******REMOVED******REMOVED*** 0.5***REMOVED******REMOVED******REMOVED*** 0.5***REMOVED***  0.522499***REMOVED***  0.544639***REMOVED***  0.566406***REMOVED***  0.587785***REMOVED***  0.608761***REMOVED******REMOVED***0.62932***REMOVED***  0.649448***REMOVED***  0.669131***REMOVED***  0.688355***REMOVED***  0.707107***REMOVED***  0.725374***REMOVED***  0.743145***REMOVED***  0.760406***REMOVED***  0.777146***REMOVED***  0.793353***REMOVED***  0.809017***REMOVED***  0.824126***REMOVED***  0.838671***REMOVED******REMOVED***0.85264***REMOVED***  0.866025***REMOVED***  0.878817***REMOVED***  0.891007***REMOVED***  0.902585***REMOVED***  0.913545***REMOVED******REMOVED***0.92388***REMOVED******REMOVED***0.93358***REMOVED***  0.942641***REMOVED***  0.951057***REMOVED******REMOVED***0.95882***REMOVED***  0.965926***REMOVED******REMOVED***0.97237***REMOVED***  0.978148***REMOVED***  0.983255***REMOVED***  0.987688***REMOVED***  0.991445***REMOVED***  0.994522***REMOVED***  0.996917***REMOVED******REMOVED***0.99863***REMOVED***  0.999657***REMOVED******REMOVED******REMOVED******REMOVED***1
***REMOVED******REMOVED******REMOVED***  1***REMOVED***  0.999657***REMOVED******REMOVED***0.99863***REMOVED***  0.996917***REMOVED***  0.994522***REMOVED***  0.991445***REMOVED***  0.987688***REMOVED***  0.983255***REMOVED***  0.978148***REMOVED******REMOVED***0.97237***REMOVED***  0.965926***REMOVED******REMOVED***0.95882***REMOVED***  0.951057***REMOVED***  0.942641***REMOVED******REMOVED***0.93358***REMOVED******REMOVED***0.92388***REMOVED***  0.913545***REMOVED***  0.902585***REMOVED***  0.891007***REMOVED***  0.878817***REMOVED***  0.866025***REMOVED***  0.866025***REMOVED***  0.878817***REMOVED***  0.891007***REMOVED***  0.902585***REMOVED***  0.913545***REMOVED******REMOVED***0.92388***REMOVED******REMOVED***0.93358***REMOVED***  0.942641***REMOVED***  0.951057***REMOVED******REMOVED***0.95882***REMOVED***  0.965926***REMOVED******REMOVED***0.97237***REMOVED***  0.978148***REMOVED***  0.983255***REMOVED***  0.987688***REMOVED***  0.991445***REMOVED***  0.994522***REMOVED***  0.996917***REMOVED******REMOVED***0.99863***REMOVED***  0.999657***REMOVED******REMOVED******REMOVED******REMOVED***1***REMOVED******REMOVED******REMOVED******REMOVED***2***REMOVED******REMOVED***1.99931***REMOVED******REMOVED***1.99726***REMOVED******REMOVED***1.99383***REMOVED******REMOVED***1.98904***REMOVED******REMOVED***1.98289***REMOVED******REMOVED***1.97538***REMOVED******REMOVED***1.96651***REMOVED******REMOVED*** 1.9563***REMOVED******REMOVED***1.94474***REMOVED******REMOVED***1.93185***REMOVED******REMOVED***1.91764***REMOVED******REMOVED***1.90211***REMOVED******REMOVED***1.88528***REMOVED******REMOVED***1.86716***REMOVED******REMOVED***1.84776***REMOVED******REMOVED***1.82709***REMOVED******REMOVED***1.80517***REMOVED******REMOVED***1.78201***REMOVED******REMOVED***1.75763***REMOVED******REMOVED***1.73205***REMOVED******REMOVED******REMOVED******REMOVED***1***REMOVED***  0.999657***REMOVED******REMOVED***0.99863***REMOVED***  0.996917***REMOVED***  0.994522***REMOVED***  0.991445***REMOVED***  0.987688***REMOVED***  0.983255***REMOVED***  0.978148***REMOVED******REMOVED***0.97237***REMOVED***  0.965926***REMOVED******REMOVED***0.95882***REMOVED***  0.951057***REMOVED***  0.942641***REMOVED******REMOVED***0.93358***REMOVED******REMOVED***0.92388***REMOVED***  0.913545***REMOVED***  0.902585***REMOVED***  0.891007***REMOVED***  0.878817***REMOVED***  0.866025***REMOVED***  0.866025***REMOVED******REMOVED***0.85264***REMOVED***  0.838671***REMOVED***  0.824126***REMOVED***  0.809017***REMOVED***  0.793353***REMOVED***  0.777146***REMOVED***  0.760406***REMOVED***  0.743145***REMOVED***  0.725374***REMOVED***  0.707107***REMOVED***  0.688355***REMOVED***  0.669131***REMOVED***  0.649448***REMOVED******REMOVED***0.62932***REMOVED***  0.608761***REMOVED***  0.587785***REMOVED***  0.566406***REMOVED***  0.544639***REMOVED***  0.522499***REMOVED******REMOVED******REMOVED*** 0.5***REMOVED***  0.477159***REMOVED******REMOVED***0.45399***REMOVED***  0.430511***REMOVED***  0.406737***REMOVED***  0.382683***REMOVED***  0.358368***REMOVED***  0.333807***REMOVED***  0.309017***REMOVED***  0.284015***REMOVED***  0.258819***REMOVED***  0.233445***REMOVED***  0.207912***REMOVED***  0.182236***REMOVED***  0.156434***REMOVED***  0.130526***REMOVED***  0.104528***REMOVED*** 0.0784591***REMOVED***  0.052336***REMOVED*** 0.0261769  4.01902e-12

Independent variable at boundary nodes (tBC) =***REMOVED******REMOVED***  0 0.523599 0.523599***REMOVED***1.0472***REMOVED***1.0472***REMOVED***1.5708***REMOVED***1.5708***REMOVED***2.0944***REMOVED***2.0944  3.14159

State vector at boundary nodes (xBC) = 
 1.40454e-28***REMOVED******REMOVED******REMOVED*** 0.5***REMOVED******REMOVED******REMOVED***-0.5 -1.42748e-12 -1.42748e-12***REMOVED******REMOVED******REMOVED******REMOVED***1 -4.25982e-12***REMOVED******REMOVED******REMOVED*** 0.5***REMOVED******REMOVED******REMOVED*** 0.5***REMOVED******REMOVED******REMOVED******REMOVED***1
***REMOVED******REMOVED******REMOVED***  1***REMOVED***  0.866025***REMOVED***  0.866025***REMOVED******REMOVED******REMOVED******REMOVED***1***REMOVED******REMOVED******REMOVED******REMOVED***2***REMOVED******REMOVED***1.73205***REMOVED******REMOVED******REMOVED******REMOVED***1***REMOVED***  0.866025***REMOVED***  0.866025  4.01902e-12

====================================================================================================

========================
Boundary Value Problem 4
========================

Boundary nodes (tBC) = 
***REMOVED******REMOVED*** 0 0.785398***REMOVED***1.5708  2.35619  3.14159

Starting state vector (oxt1) = 
-0.1  0.7  0.1 -0.7
 1.1  0.7 -1.1 -0.7

Initiating the BVP solver...

Boundary nodes correspond to the below columns: 
  0  25  50  75 100

k = 0... kG = 0.256924... kGPrev = inf alpha = 1 Increasing alpha...
k = 1... kG = 0.0456267... kGPrev = 0.256924 alpha = 1
k = 2... kG = 0.0111092... kGPrev = 0.0456267 alpha = 1
k = 3... kG = 0.00239144... kGPrev = 0.0111092 alpha = 1
k = 4... kG = 0.000231811... kGPrev = 0.00239144 alpha = 1
k = 5... kG = 2.63355e-06... kGPrev = 0.000231811 alpha = 1
k = 6... kG = 3.41465e-10... kGPrev = 2.63355e-06 alpha = 1 Increasing alpha...
Ran 7 iteration(s).

========
Solution
========
Independent variable (t) = 
***REMOVED******REMOVED***  0 0.0314159 0.0628319 0.0942478  0.125664***REMOVED***0.15708  0.188496  0.219911  0.251327  0.282743  0.314159  0.345575  0.376991  0.408407  0.439823  0.471239  0.502655  0.534071  0.565487  0.596903  0.628319  0.659734***REMOVED***0.69115  0.722566  0.753982  0.785398  0.785398  0.816814***REMOVED***0.84823  0.879646  0.911062  0.942478  0.973894***REMOVED***1.00531***REMOVED***1.03673***REMOVED***1.06814***REMOVED***1.09956***REMOVED***1.13097***REMOVED***1.16239***REMOVED***1.19381***REMOVED***1.22522***REMOVED***1.25664***REMOVED***1.28805***REMOVED***1.31947***REMOVED***1.35088***REMOVED*** 1.3823***REMOVED***1.41372***REMOVED***1.44513***REMOVED***1.47655***REMOVED***1.50796***REMOVED***1.53938***REMOVED*** 1.5708***REMOVED*** 1.5708***REMOVED***1.60221***REMOVED***1.63363***REMOVED***1.66504***REMOVED***1.69646***REMOVED***1.72788***REMOVED***1.75929***REMOVED***1.79071***REMOVED***1.82212***REMOVED***1.85354***REMOVED***1.88496***REMOVED***1.91637***REMOVED***1.94779***REMOVED*** 1.9792***REMOVED***2.01062***REMOVED***2.04204***REMOVED***2.07345***REMOVED***2.10487***REMOVED***2.13628***REMOVED*** 2.1677***REMOVED***2.19911***REMOVED***2.23053***REMOVED***2.26195***REMOVED***2.29336***REMOVED***2.32478***REMOVED***2.35619***REMOVED***2.35619***REMOVED***2.38761***REMOVED***2.41903***REMOVED***2.45044***REMOVED***2.48186***REMOVED***2.51327***REMOVED***2.54469***REMOVED***2.57611***REMOVED***2.60752***REMOVED***2.63894***REMOVED***2.67035***REMOVED***2.70177***REMOVED***2.73319***REMOVED*** 2.7646***REMOVED***2.79602***REMOVED***2.82743***REMOVED***2.85885***REMOVED***2.89027***REMOVED***2.92168***REMOVED*** 2.9531***REMOVED***2.98451***REMOVED***3.01593***REMOVED***3.04734***REMOVED***3.07876***REMOVED***3.11018***REMOVED***3.14159

State vector (x) = 
 2.16514e-28***REMOVED*** 0.0314108***REMOVED*** 0.0627905***REMOVED*** 0.0941083***REMOVED***  0.125333***REMOVED***  0.156434***REMOVED***  0.187381***REMOVED***  0.218143***REMOVED******REMOVED***0.24869***REMOVED***  0.278991***REMOVED***  0.309017***REMOVED***  0.338738***REMOVED***  0.368125***REMOVED***  0.397148***REMOVED***  0.425779***REMOVED******REMOVED***0.45399***REMOVED***  0.481754***REMOVED***  0.509041***REMOVED***  0.535827***REMOVED***  0.562083***REMOVED***  0.587785***REMOVED***  0.612907***REMOVED***  0.637424***REMOVED***  0.661312***REMOVED***  0.684547***REMOVED***  0.707107***REMOVED***  0.707107***REMOVED***  0.728969***REMOVED***  0.750111***REMOVED***  0.770513***REMOVED***  0.790155***REMOVED***  0.809017***REMOVED***  0.827081***REMOVED***  0.844328***REMOVED***  0.860742***REMOVED***  0.876307***REMOVED***  0.891007***REMOVED***  0.904827***REMOVED***  0.917755***REMOVED***  0.929776***REMOVED***  0.940881***REMOVED***  0.951057***REMOVED***  0.960294***REMOVED***  0.968583***REMOVED***  0.975917***REMOVED***  0.982287***REMOVED***  0.987688***REMOVED***  0.992115***REMOVED***  0.995562***REMOVED***  0.998027***REMOVED***  0.999507***REMOVED******REMOVED******REMOVED******REMOVED***1***REMOVED*** -0.257684***REMOVED***  -0.26603***REMOVED*** -0.274114***REMOVED*** -0.281927***REMOVED*** -0.289463***REMOVED*** -0.296712***REMOVED*** -0.303669***REMOVED*** -0.310326***REMOVED*** -0.316676***REMOVED*** -0.322715***REMOVED*** -0.328434***REMOVED***  -0.33383***REMOVED*** -0.338896***REMOVED*** -0.343628***REMOVED*** -0.348021***REMOVED***  -0.35207***REMOVED*** -0.355771***REMOVED*** -0.359122***REMOVED*** -0.362118***REMOVED*** -0.364757***REMOVED*** -0.367036***REMOVED*** -0.368952***REMOVED*** -0.370505***REMOVED*** -0.371692***REMOVED*** -0.372512***REMOVED*** -0.372964***REMOVED*** -0.372964***REMOVED*** -0.373049***REMOVED*** -0.372765***REMOVED*** -0.372113***REMOVED*** -0.371094***REMOVED*** -0.369709***REMOVED*** -0.367959***REMOVED*** -0.365846***REMOVED*** -0.363372***REMOVED*** -0.360539***REMOVED*** -0.357351***REMOVED*** -0.353809***REMOVED*** -0.349919***REMOVED*** -0.345683***REMOVED*** -0.341106***REMOVED*** -0.336193***REMOVED*** -0.330948***REMOVED*** -0.325376***REMOVED*** -0.319483***REMOVED*** -0.313274***REMOVED*** -0.306757***REMOVED*** -0.299937***REMOVED*** -0.292821***REMOVED*** -0.285415***REMOVED*** -0.277729***REMOVED*** -0.269768
***REMOVED******REMOVED******REMOVED***  1***REMOVED***  0.999507***REMOVED***  0.998027***REMOVED***  0.995562***REMOVED***  0.992115***REMOVED***  0.987688***REMOVED***  0.982287***REMOVED***  0.975917***REMOVED***  0.968583***REMOVED***  0.960294***REMOVED***  0.951057***REMOVED***  0.940881***REMOVED***  0.929776***REMOVED***  0.917755***REMOVED***  0.904827***REMOVED***  0.891007***REMOVED***  0.876307***REMOVED***  0.860742***REMOVED***  0.844328***REMOVED***  0.827081***REMOVED***  0.809017***REMOVED***  0.790155***REMOVED***  0.770513***REMOVED***  0.750111***REMOVED***  0.728969***REMOVED***  0.707107***REMOVED***  0.707107***REMOVED***  0.684547***REMOVED***  0.661312***REMOVED***  0.637424***REMOVED***  0.612907***REMOVED***  0.587785***REMOVED***  0.562083***REMOVED***  0.535827***REMOVED***  0.509041***REMOVED***  0.481754***REMOVED******REMOVED***0.45399***REMOVED***  0.425779***REMOVED***  0.397148***REMOVED***  0.368125***REMOVED***  0.338738***REMOVED***  0.309017***REMOVED***  0.278991***REMOVED******REMOVED***0.24869***REMOVED***  0.218143***REMOVED***  0.187381***REMOVED***  0.156434***REMOVED***  0.125333***REMOVED*** 0.0941083***REMOVED*** 0.0627905***REMOVED*** 0.0314108 -7.18934e-13***REMOVED*** -0.269768***REMOVED***  -0.26154***REMOVED*** -0.253055***REMOVED***  -0.24432***REMOVED*** -0.235344***REMOVED*** -0.226136***REMOVED*** -0.216704***REMOVED*** -0.207059***REMOVED*** -0.197209***REMOVED*** -0.187165***REMOVED*** -0.176936***REMOVED*** -0.166532***REMOVED*** -0.155964***REMOVED*** -0.145242***REMOVED*** -0.134377***REMOVED*** -0.123379***REMOVED*** -0.112259***REMOVED*** -0.101029***REMOVED***-0.0896985***REMOVED***-0.0782798***REMOVED***-0.0667839***REMOVED***-0.0552221***REMOVED***-0.0436058***REMOVED***-0.0319464***REMOVED***-0.0202555  -0.00854464  -0.00854464***REMOVED***0.00317467***REMOVED*** 0.0148908***REMOVED*** 0.0265923***REMOVED*** 0.0382676***REMOVED***  0.049905***REMOVED*** 0.0614933***REMOVED*** 0.0730208***REMOVED*** 0.0844763***REMOVED*** 0.0958484***REMOVED***  0.107126***REMOVED***  0.118298***REMOVED***  0.129353***REMOVED******REMOVED***0.14028***REMOVED***  0.151069***REMOVED***  0.161709***REMOVED***  0.172189***REMOVED******REMOVED*** 0.1825***REMOVED******REMOVED***0.19263***REMOVED******REMOVED***0.20257***REMOVED******REMOVED***0.21231***REMOVED***  0.221841***REMOVED***  0.231153***REMOVED***  0.240236***REMOVED***  0.249083***REMOVED***  0.257684

Independent variable at boundary nodes (tBC) =***REMOVED******REMOVED***  0 0.785398 0.785398***REMOVED***1.5708***REMOVED***1.5708  2.35619  2.35619  3.14159

State vector at boundary nodes (xBC) = 
 2.16514e-28***REMOVED***  0.707107***REMOVED***  0.707107***REMOVED******REMOVED******REMOVED******REMOVED***1***REMOVED*** -0.257684***REMOVED*** -0.372964***REMOVED*** -0.372964***REMOVED*** -0.269768
***REMOVED******REMOVED******REMOVED***  1***REMOVED***  0.707107***REMOVED***  0.707107 -7.18934e-13***REMOVED*** -0.269768  -0.00854464  -0.00854464***REMOVED***  0.257684

====================================================================================================

Program ended...

[shivanandvp@computer nlmpbvp]$
```

## Source - test_nlmpBVP.cpp:
```cpp
// ========================================
// Author: shivanandvp 
// Email : shivanandvp.oss@gmail.com
// ========================================
// Copyright shivanandvp (shivanandvp.oss@gmail.com)

// ==========
// References
// ==========
// [1] Ojika, T., and Y. Kasue. "Initial-value adjusting method for the solution of nonlinear multipoint boundary-value problems."
//***REMOVED***  Journal of Mathematical Analysis and Applications 69.2 (1979): 359-371.
// [2] Welsh, Wayne, and Takeo Ojika. "Multipoint boundary value problems with discontinuities I. Algorithms and applications."
//***REMOVED***  Journal of Computational and Applied Mathematics 6.2 (1980): 133-143.

// ============================
// Credits for example problems
// ============================
/* Boundary Value Problem 1 */
// "Using Initial Guess to Indicate Desired Solution"
// “bvp4c.” Mathworks Documentation, Mathworks, www.mathworks.com/help/matlab/ref/bvp4c.html

/* Boundary value problem 2 */
// Example 2
// Tr, Ramesh. (2017). A novel method for solving multipoint boundary value problems.
// Global Journal of Pure and Applied Mathematics. 13. 850-857. 

/* Boundary value problem 3 */
// Example 1
// Welsh, Wayne, and Takeo Ojika.
// "Multipoint boundary value problems with discontinuities I. Algorithms and applications."
// Journal of Computational and Applied Mathematics 6.2 (1980): 133-143.

/* Boundary value problem 4 */
// Example 3
// Welsh, Wayne, and Takeo Ojika.
// "Multipoint boundary value problems with discontinuities I. Algorithms and applications."
// Journal of Computational and Applied Mathematics 6.2 (1980): 133-143.
// ============================


// ===============================
// Includes and global definitions
// ===============================
#include <iostream>***REMOVED******REMOVED******REMOVED******REMOVED*** // For the cout statements
#include "nlmpbvp.hpp"***REMOVED******REMOVED******REMOVED*** // For the boundary value problem solver function declarations
#include <Eigen/MPRealSupport>  // For arbitrary precision computation
using namespace std;***REMOVED******REMOVED******REMOVED******REMOVED***// For cout
using namespace Eigen;***REMOVED******REMOVED******REMOVED*** // For matrix and vector data types and operations
using namespace mpfr;***REMOVED******REMOVED******REMOVED***  // For arbitrary precision computation
// ===============================

// ================
// Functions dxBydt
// ================
// dxBydt = a function that defines the derivative of a state vector x at t -- (nx1)

/* Boundary Value Problem 1 */
VectorXm<mpreal> dxBydt_BVP1(mpreal t, VectorXm<mpreal> x){ 
***REMOVED*** VectorXm<mpreal> dxdt(2);
***REMOVED*** dxdt(0) = x(1);
***REMOVED*** dxdt(1) = -fabs(x(0));
***REMOVED*** return dxdt;
}

/* Boundary Value Problem 2 */
VectorXm<mpreal> dxBydt_BVP2(mpreal t, VectorXm<mpreal> x){ 
***REMOVED*** VectorXm<mpreal> dxdt(3);
***REMOVED*** dxdt(0) = x(1);
***REMOVED*** dxdt(1) = x(2);
***REMOVED*** dxdt(2) = 25*x(1) - 1;
***REMOVED*** return dxdt;
}

/* Boundary Value Problem 3 */
VectorXm<mpreal> dxBydt_BVP3(mpreal t, VectorXm<mpreal> x){ 
***REMOVED*** VectorXm<mpreal> dxdt(2);
***REMOVED*** dxdt(0) =  x(1);
***REMOVED*** dxdt(1) = -x(0);
***REMOVED*** return dxdt;
}

/* Boundary Value Problem 4 */
VectorXm<mpreal> dxBydt_BVP4(mpreal t, VectorXm<mpreal> x){ 
***REMOVED*** VectorXm<mpreal> dxdt(2);
***REMOVED*** dxdt(0) =  x(1);
***REMOVED*** dxdt(1) = -x(0);
***REMOVED*** return dxdt;
}
// ================

// ====================
// Functions BCResidues
// ====================
// BCResidues = a function that defines the boundary condition residues at nodal state vectors xBC -- (nx1) 

/* Boundary Value Problem 1 */
VectorXm<mpreal> BCResidues_BVP1(MatrixXm<mpreal> xBC){
***REMOVED*** VectorXm<mpreal> residues(2);
***REMOVED*** residues(0) = xBC(0,0) - 0;
***REMOVED*** residues(1) = xBC(0,1) + 2;
***REMOVED*** return residues;
}

/* Boundary Value Problem 2 */
VectorXm<mpreal> BCResidues_BVP2(MatrixXm<mpreal> xBC){
***REMOVED*** VectorXm<mpreal> residues(3);
***REMOVED*** residues(0) = xBC(1,0) - 0;
***REMOVED*** residues(1) = xBC(1,2) - 0;
***REMOVED*** residues(2) = xBC(0,1) - 0;
***REMOVED*** return residues;
}


// BCResidues***REMOVED***  = a function that defines the boundary condition residues at the nodal state vectors... -- (n(m-1)x1)
//***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***... on the left and right side of integration intervals, xBCL and xBCR

/* Boundary Value Problem 3 */
VectorXm<mpreal> BCResidues_BVP3(MatrixXm<mpreal> xBCL, MatrixXm<mpreal> xBCR){
***REMOVED*** VectorXm<mpreal> residues(10);
***REMOVED*** residues(0) = xBCL(0,0)-0;
***REMOVED*** residues(1) = xBCL(1,0)-1;
***REMOVED*** residues(2) = xBCR(0,0) - xBCL(0,1) - 1;
***REMOVED*** residues(3) = xBCR(1,0) - xBCL(1,1) - 0;
***REMOVED*** residues(4) = xBCR(1,1) - xBCL(1,2) + 1;
***REMOVED*** residues(5) = xBCR(0,1) - xBCL(0,2) + 0;
***REMOVED*** residues(6) = xBCR(0,2) - xBCL(0,3) - 1;
***REMOVED*** residues(7) = xBCR(1,2) - xBCL(1,3) - (sqrt(3)-1);
***REMOVED*** residues(8) = xBCR(0,3) - xBCL(0,4) - 0;
***REMOVED*** residues(9) = xBCR(1,3) - xBCL(1,4) - 0;
***REMOVED*** return residues;
}

/* Boundary Value Problem 4 */
VectorXm<mpreal> BCResidues_BVP4(MatrixXm<mpreal> xBCL, MatrixXm<mpreal> xBCR){
***REMOVED*** VectorXm<mpreal> residues(8);
***REMOVED*** residues(0) = xBCL(0,0)-0;
***REMOVED*** residues(1) = xBCL(1,0)-1;
***REMOVED*** residues(2) = pow(xBCL(0,1),2)*pow(xBCL(1,1),3) - exp(xBCR(1,1))*pow(xBCR(1,0),2)*pow(xBCL(1,3),2) + pow(xBCR(1,3),2)*pow(xBCR(0,2),2) - 0.1859767072;
***REMOVED*** residues(3) = pow(xBCR(0,1),2)*pow(xBCR(0,0),3)*pow(xBCL(0,3),2) + pow(xBCL(1,2),2)*exp(xBCR(1,2)) + pow(xBCL(0,2),2)*pow(xBCR(0,3),2) - 0.1261677772;
***REMOVED*** residues(4) = xBCR(0,0) - xBCL(0,1) + 0;
***REMOVED*** residues(5) = xBCR(1,0) - xBCL(1,1) + 0;
***REMOVED*** residues(6) = xBCR(0,2) - xBCL(0,3) + 0;
***REMOVED*** residues(7) = xBCR(1,2) - xBCL(1,3) + 0;
***REMOVED*** return residues;
}


// =================
// The main function
// =================
// This is where the program execution begins
int main(
***REMOVED*** int argc,***REMOVED***// argc = the number of arguments passed to the program
***REMOVED*** char **argv // argv = an array of strings that are passed to  the program 
***REMOVED*** ){

***REMOVED*** mpreal::set_default_prec(64); // Set the number of digits of precision you want for computations

***REMOVED*** cout<<endl;
***REMOVED*** cout<<"============================================================================="<<endl;
***REMOVED*** cout<<"Test: Non-linear multipoint boundary value problem solver (nlmpBVP, nlmpBVP2)"<<endl;
***REMOVED*** cout<<"============================================================================="<<endl;
***REMOVED*** cout<<"Copyright shivanandvp (shivanandvp.oss@gmail.com)"<<endl;

***REMOVED*** // Variable declarations***REMOVED***

***REMOVED*** /* Boundary Value Problem 1 */
***REMOVED*** RowVectorXm<mpreal> tBC_BVP1(2);***REMOVED******REMOVED******REMOVED******REMOVED***// t_BC***REMOVED******REMOVED******REMOVED***  = row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED*** VectorXm<mpreal>***REMOVED***oxt1_BVP1(2);***REMOVED******REMOVED******REMOVED******REMOVED***// oxt1***REMOVED******REMOVED******REMOVED***  = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)

***REMOVED*** /* Boundary Value Problem 2 */
***REMOVED*** RowVectorXm<mpreal> tBC_BVP2(3);***REMOVED******REMOVED******REMOVED******REMOVED***// t_BC***REMOVED******REMOVED******REMOVED***  = row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED*** VectorXm<mpreal>***REMOVED***oxt1_BVP2(3);***REMOVED******REMOVED******REMOVED******REMOVED***// oxt1***REMOVED******REMOVED******REMOVED***  = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***-- (nx1)

***REMOVED*** /* Boundary Value Problem 3 */
***REMOVED*** RowVectorXm<mpreal> tBC_BVP3(6);***REMOVED******REMOVED******REMOVED******REMOVED***// t_BC***REMOVED******REMOVED******REMOVED***  = row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED******REMOVED***  -- (1xm)  
***REMOVED*** MatrixXm<mpreal> oxt1_BVP3(2,5);***REMOVED******REMOVED******REMOVED******REMOVED***// oxt1***REMOVED******REMOVED******REMOVED***  = matrix of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx(m-1))

***REMOVED*** /* Boundary Value Problem 4 */
***REMOVED*** RowVectorXm<mpreal> tBC_BVP4(5);***REMOVED******REMOVED******REMOVED******REMOVED***// t_BC***REMOVED******REMOVED******REMOVED***  = row vector of values at which the boundary conditions are specified***REMOVED******REMOVED******REMOVED******REMOVED***  -- (1xm)
***REMOVED*** MatrixXm<mpreal> oxt1_BVP4(2,4);***REMOVED******REMOVED******REMOVED******REMOVED***// oxt1***REMOVED******REMOVED******REMOVED***  = matrix of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx(m-1))

***REMOVED*** BVPSolution<mpreal> bvpSolution_BVP1;***REMOVED******REMOVED*** // bvpSolution***REMOVED*** = the structure in which the solutions of the boundary value problem will be saved
***REMOVED*** IVAMParameters<mpreal> ivamParameters_BVP1; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

***REMOVED*** BVPSolution<mpreal> bvpSolution_BVP2;***REMOVED******REMOVED*** // bvpSolution***REMOVED*** = the structure in which the solutions of the boundary value problem will be saved
***REMOVED*** IVAMParameters<mpreal> ivamParameters_BVP2; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

***REMOVED*** BVPSolution<mpreal> bvpSolution_BVP3;***REMOVED******REMOVED*** // bvpSolution***REMOVED*** = the structure in which the solutions of the boundary value problem will be saved
***REMOVED*** IVAMParameters<mpreal> ivamParameters_BVP3; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

***REMOVED*** BVPSolution<mpreal> bvpSolution_BVP4;***REMOVED******REMOVED*** // bvpSolution***REMOVED*** = the structure in which the solutions of the boundary value problem will be saved
***REMOVED*** IVAMParameters<mpreal> ivamParameters_BVP4; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

***REMOVED*** // Variable definitions

***REMOVED*** /* Boundary Value Problem 1 */
***REMOVED*** tBC_BVP1  << 0.0, 4.0;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
***REMOVED*** oxt1_BVP1 <<  1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// oxt1 = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx1) 
***REMOVED******REMOVED******REMOVED******REMOVED*** 0;  
***REMOVED*** // tBC_BVP1  << 0.0, 4.0;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
***REMOVED*** // oxt1_BVP1 << -1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED***// oxt1 = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx1) 
***REMOVED*** //***REMOVED******REMOVED******REMOVED*** 0;

***REMOVED*** /* Boundary Value Problem 2 */
***REMOVED*** tBC_BVP2  << 0.0, 0.5, 1.0;***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
***REMOVED*** oxt1_BVP2 << 1,***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** // oxt1 = column vector of the guessed initial state***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx1)
***REMOVED******REMOVED******REMOVED******REMOVED***1,
***REMOVED******REMOVED******REMOVED******REMOVED***1;

***REMOVED*** /* Boundary Value Problem 3 */
***REMOVED*** // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
***REMOVED*** tBC_BVP3  << 0.0, mpfr::const_pi()/6, mpfr::const_pi()/3, mpfr::const_pi()/2, 2*mpfr::const_pi()/3, mpfr::const_pi();
***REMOVED*** // oxt1 = a matrix of the guessed initial state on the left side of each integration interval***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx(m-1))
***REMOVED*** oxt1_BVP3 <<  0.1, 0.1, 0.4, 0.8, 0.9,
***REMOVED******REMOVED******REMOVED******REMOVED***-0.6, 0.1, 0.9, 2.1, 0.8;

***REMOVED*** /* Boundary Value Problem 4 */
***REMOVED*** // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
***REMOVED*** tBC_BVP4  << 0.0, mpfr::const_pi()/4, mpfr::const_pi()/2, 3*mpfr::const_pi()/4, mpfr::const_pi();
***REMOVED*** // oxt1 = a matrix of the guessed initial state on the left side of each integration interval***REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED******REMOVED*** -- (nx(m-1))
***REMOVED*** oxt1_BVP4 << -0.1, 0.7, 0.1,-0.7,
***REMOVED******REMOVED******REMOVED******REMOVED*** 1.1, 0.7,-1.1,-0.7;

***REMOVED*** // Assign the parameters for IVAM

***REMOVED*** ivamParameters_BVP1.EPSILON***REMOVED*** = 1e-10; // EPSILON***REMOVED*** = the state perturbation parameter to probe the differential equation system with
***REMOVED*** ivamParameters_BVP1.ALPHA***REMOVED******REMOVED***= 1.0;***REMOVED***// ALPHA***REMOVED******REMOVED***= the relaxation factor to scale the adjustment to the initial condition
***REMOVED*** ivamParameters_BVP1.SIGMA***REMOVED******REMOVED***= 1e-14; // SIGMA***REMOVED******REMOVED***= the tolerance for error outside which the solver needs to  iterate further. 
***REMOVED*** ivamParameters_BVP1.BETA***REMOVED******REMOVED*** = 1e-3;  // BETA***REMOVED******REMOVED*** = the deflation factor
***REMOVED*** ivamParameters_BVP1.printDebug = true;  // printDebug = specify whether debug messages should be output to the console

***REMOVED*** ivamParameters_BVP2.EPSILON***REMOVED*** = 1e-10; // EPSILON***REMOVED*** = the state perturbation parameter to probe the differential equation system with
***REMOVED*** ivamParameters_BVP2.ALPHA***REMOVED******REMOVED***= 1.0;***REMOVED***// ALPHA***REMOVED******REMOVED***= the relaxation factor to scale the adjustment to the initial condition
***REMOVED*** ivamParameters_BVP2.SIGMA***REMOVED******REMOVED***= 1e-14; // SIGMA***REMOVED******REMOVED***= the tolerance for error outside which the solver needs to  iterate further. 
***REMOVED*** ivamParameters_BVP2.BETA***REMOVED******REMOVED*** = 1e-3;  // BETA***REMOVED******REMOVED*** = the deflation factor
***REMOVED*** ivamParameters_BVP2.printDebug = true;  // printDebug = specify whether debug messages should be output to the console

***REMOVED*** ivamParameters_BVP3.EPSILON***REMOVED*** = 1e-10; // EPSILON***REMOVED*** = the state perturbation parameter to probe the differential equation system with
***REMOVED*** ivamParameters_BVP3.ALPHA***REMOVED******REMOVED***= 1.0;***REMOVED***// ALPHA***REMOVED******REMOVED***= the relaxation factor to scale the adjustment to the initial condition
***REMOVED*** ivamParameters_BVP3.SIGMA***REMOVED******REMOVED***= 1e-14; // SIGMA***REMOVED******REMOVED***= the tolerance for error outside which the solver needs to  iterate further. 
***REMOVED*** ivamParameters_BVP3.BETA***REMOVED******REMOVED*** = 1e-3;  // BETA***REMOVED******REMOVED*** = the deflation factor
***REMOVED*** ivamParameters_BVP3.printDebug = true;  // printDebug = specify whether debug messages should be output to the console

***REMOVED*** ivamParameters_BVP4.EPSILON***REMOVED*** = 1e-10; // EPSILON***REMOVED*** = the state perturbation parameter to probe the differential equation system with
***REMOVED*** ivamParameters_BVP4.ALPHA***REMOVED******REMOVED***= 1.0;***REMOVED***// ALPHA***REMOVED******REMOVED***= the relaxation factor to scale the adjustment to the initial condition
***REMOVED*** ivamParameters_BVP4.SIGMA***REMOVED******REMOVED***= 1e-14; // SIGMA***REMOVED******REMOVED***= the tolerance for error outside which the solver needs to  iterate further. 
***REMOVED*** ivamParameters_BVP4.BETA***REMOVED******REMOVED*** = 1e-3;  // BETA***REMOVED******REMOVED*** = the deflation factor
***REMOVED*** ivamParameters_BVP4.printDebug = true;  // printDebug = specify whether debug messages should be output to the console

***REMOVED*** /* Boundary Value Problem 1 */
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========================"<<endl;
***REMOVED*** cout<<"Boundary Value Problem 1"<<endl;
***REMOVED*** cout<<"========================"<<endl;
***REMOVED*** cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC_BVP1<<endl;
***REMOVED*** cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1_BVP1<<endl;
***REMOVED*** cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
***REMOVED*** bvpSolution_BVP1 = nlmpBVP<mpreal>(2, 2, 101, tBC_BVP1, oxt1_BVP1, dxBydt_BVP1, BCResidues_BVP1, ivamParameters_BVP1);
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Solution"<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Independent variable (t) = ";
***REMOVED*** cout<<endl<<bvpSolution_BVP1.t<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector (x) = "<<endl;
***REMOVED*** cout<<bvpSolution_BVP1.x<<endl;
***REMOVED*** cout<<endl;***REMOVED*** 
***REMOVED*** cout<<"Independent variable at boundary nodes (tBC) = ";
***REMOVED*** cout<<bvpSolution_BVP1.tBC<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector at boundary nodes (xBC) = "<<endl;
***REMOVED*** cout<<bvpSolution_BVP1.xBC<<endl;
***REMOVED*** cout<<endl; 
***REMOVED*** cout<<"===================================================================================================="<<endl; 

***REMOVED***  /* Boundary Value Problem 2 */
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========================"<<endl;
***REMOVED*** cout<<"Boundary Value Problem 2"<<endl;
***REMOVED*** cout<<"========================"<<endl;
***REMOVED*** cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC_BVP2<<endl;
***REMOVED*** cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1_BVP2<<endl;
***REMOVED*** cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
***REMOVED*** bvpSolution_BVP2 = nlmpBVP<mpreal>(3, 3, 101, tBC_BVP2, oxt1_BVP2, dxBydt_BVP2, BCResidues_BVP2, ivamParameters_BVP2);
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Solution"<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Independent variable (t) = ";
***REMOVED*** cout<<endl<<bvpSolution_BVP2.t<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector (x) = "<<endl;
***REMOVED*** cout<<bvpSolution_BVP2.x<<endl;
***REMOVED*** cout<<endl;***REMOVED*** 
***REMOVED*** cout<<"Independent variable at boundary nodes (tBC) = ";
***REMOVED*** cout<<bvpSolution_BVP2.tBC<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector at boundary nodes (xBC) = "<<endl;
***REMOVED*** cout<<bvpSolution_BVP2.xBC<<endl;
***REMOVED*** cout<<endl; 
***REMOVED*** cout<<"===================================================================================================="<<endl; 

***REMOVED*** /* Boundary Value Problem 3 */
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========================"<<endl;
***REMOVED*** cout<<"Boundary Value Problem 3"<<endl;
***REMOVED*** cout<<"========================"<<endl;
***REMOVED*** cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC_BVP3<<endl;
***REMOVED*** cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1_BVP3<<endl;
***REMOVED*** cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
***REMOVED*** bvpSolution_BVP3 = nlmpBVP2<mpreal>(2, 6, 12*10+1, tBC_BVP3, oxt1_BVP3, dxBydt_BVP3, BCResidues_BVP3, ivamParameters_BVP3);
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Solution"<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Independent variable (t) = ";
***REMOVED*** cout<<endl<<bvpSolution_BVP3.t<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector (x) = "<<endl;
***REMOVED*** cout<<bvpSolution_BVP3.x<<endl;
***REMOVED*** cout<<endl;***REMOVED*** 
***REMOVED*** cout<<"Independent variable at boundary nodes (tBC) = ";
***REMOVED*** cout<<bvpSolution_BVP3.tBC<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector at boundary nodes (xBC) = "<<endl;
***REMOVED*** cout<<bvpSolution_BVP3.xBC<<endl;
***REMOVED*** cout<<endl; 
***REMOVED*** cout<<"===================================================================================================="<<endl; 

***REMOVED*** /* Boundary Value Problem 4 */
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========================"<<endl;
***REMOVED*** cout<<"Boundary Value Problem 4"<<endl;
***REMOVED*** cout<<"========================"<<endl;
***REMOVED*** cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC_BVP4<<endl;
***REMOVED*** cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1_BVP4<<endl;
***REMOVED*** cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
***REMOVED*** bvpSolution_BVP4 = nlmpBVP2<mpreal>(2, 5, 10*10+1, tBC_BVP4, oxt1_BVP4, dxBydt_BVP4, BCResidues_BVP4, ivamParameters_BVP4);
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Solution"<<endl;
***REMOVED*** cout<<"========"<<endl;
***REMOVED*** cout<<"Independent variable (t) = ";
***REMOVED*** cout<<endl<<bvpSolution_BVP4.t<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector (x) = "<<endl;
***REMOVED*** cout<<bvpSolution_BVP4.x<<endl;
***REMOVED*** cout<<endl;***REMOVED*** 
***REMOVED*** cout<<"Independent variable at boundary nodes (tBC) = ";
***REMOVED*** cout<<bvpSolution_BVP4.tBC<<endl;
***REMOVED*** cout<<endl;
***REMOVED*** cout<<"State vector at boundary nodes (xBC) = "<<endl;
***REMOVED*** cout<<bvpSolution_BVP4.xBC<<endl;
***REMOVED*** cout<<endl; 
***REMOVED*** cout<<"===================================================================================================="<<endl; 

***REMOVED*** cout<<endl<<"Program ended..."<<endl<<endl;

***REMOVED*** return 0;
}
// =================
```

