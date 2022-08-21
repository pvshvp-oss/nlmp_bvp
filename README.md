# nlmp_bvp

**Nonlinear Multipoint Boundary Value Problem Solver** in `C++` based on:  
1. Ojika, T., and Y. Kasue. "Initial-value adjusting method for the solution of nonlinear multipoint boundary-value problems."  
Journal of Mathematical Analysis and Applications 69.2 (1979): 359-371.  
2. Welsh, Wayne, and Takeo Ojika. "Multipoint boundary value problems with discontinuities I. Algorithms and applications."  
Journal of Computational and Applied Mathematics 6.2 (1980): 133-143.  
  
This is a tool consisting of **multiple header files** written by [shivanandvp](https://www.github.com/shivanandvp) in response to the need for a suitable C++ alternative to bvp4c in Mathworks MATLAB.  
You can use the `nlmpbvp.hpp` header file from this project to be able to call the functions `nlmpBVP` and `nlmpBVP2`.  
However, this header file depends on other files too, so if you want to use this tool in your own projects, please either modify `test_nlmpBVP.cpp` or refer to `CMakeLists.txt` to find dependencies,

>>>
Note: The file `test_nlmpBVP.cpp` is only a sample file where you call the functions that are actually defined in the header files found in `include`
>>>

## Dependencies

Eigen, MPFR, and MPReal. In Arch Linux, install the packages `eigen`, `mpfr`, and `mpreal-git`.

## Output of sample user-side program - test_nlmpBVP.cpp:
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
   0 0.04 0.08 0.12 0.16  0.2 0.24 0.28 0.32 0.36  0.4 0.44 0.48 0.52 0.56  0.6 0.64 0.68 0.72 0.76  0.8 0.84 0.88 0.92 0.96    1 1.04 1.08 1.12 1.16  1.2 1.24 1.28 1.32 1.36  1.4 1.44 1.48 1.52 1.56  1.6 1.64 1.68 1.72 1.76  1.8 1.84 1.88 1.92 1.96    2 2.04 2.08 2.12 2.16  2.2 2.24 2.28 2.32 2.36  2.4 2.44 2.48 2.52 2.56  2.6 2.64 2.68 2.72 2.76  2.8 2.84 2.88 2.92 2.96    3 3.04 3.08 3.12 3.16  3.2 3.24 3.28 3.32 3.36  3.4 3.44 3.48 3.52 3.56  3.6 3.64 3.68 3.72 3.76  3.8 3.84 3.88 3.92 3.96    4

State vector (x) = 
          0 -0.00293227 -0.00586923 -0.00881558   -0.011776  -0.0147553  -0.0177583  -0.0207896  -0.0238542   -0.026957  -0.0301029  -0.0332969  -0.0365443  -0.0398501  -0.0432198  -0.0466585   -0.050172  -0.0537657  -0.0574454  -0.0612171  -0.0650867  -0.0690605  -0.0731448  -0.0773462  -0.0816713  -0.0861271  -0.0907208  -0.0954596   -0.100351   -0.105403   -0.110624   -0.116022   -0.121605   -0.127384   -0.133365   -0.139561   -0.145979   -0.152632   -0.159528    -0.16668   -0.174099   -0.181796   -0.189784   -0.198075   -0.206684   -0.215624   -0.224908   -0.234552   -0.244572   -0.254983   -0.265802   -0.277047   -0.288734   -0.300884   -0.313516   -0.326648   -0.340304   -0.354504   -0.369272    -0.38463   -0.400604   -0.417219   -0.434502    -0.45248   -0.471182   -0.490638   -0.510879   -0.531937   -0.553847   -0.576643   -0.600362   -0.625042   -0.650722   -0.677443   -0.705248   -0.734181    -0.76429   -0.795621   -0.828226   -0.862156   -0.897466   -0.934212   -0.972452    -1.01225    -1.05367    -1.09677    -1.14163    -1.18831     -1.2369    -1.28746    -1.34009    -1.39486    -1.45186    -1.51119    -1.57293    -1.63719    -1.70407    -1.77368    -1.84613    -1.92153          -2
 -0.0732871  -0.0733458  -0.0735218  -0.0738154  -0.0742272  -0.0747578   -0.075408  -0.0761788  -0.0770716  -0.0780877  -0.0792287  -0.0804965  -0.0818932  -0.0834209   -0.085082  -0.0868794  -0.0888157  -0.0908942  -0.0931181   -0.095491  -0.0980168   -0.100699   -0.103543   -0.106553   -0.109732   -0.113088   -0.116624   -0.120348   -0.124263   -0.128378   -0.132698    -0.13723   -0.141982   -0.146961   -0.152175   -0.157633   -0.163343   -0.169315   -0.175557    -0.18208   -0.188895   -0.196012   -0.203443   -0.211199   -0.219293   -0.227738   -0.236547   -0.245735   -0.255316   -0.265306   -0.275721   -0.286576    -0.29789   -0.309681   -0.321967   -0.334769   -0.348106   -0.362001   -0.376474    -0.39155   -0.407253   -0.423607   -0.440639   -0.458376   -0.476847   -0.496081   -0.516109   -0.536962   -0.558675   -0.581282   -0.604819   -0.629324   -0.654836   -0.681395   -0.709045    -0.73783   -0.767796    -0.79899   -0.831462   -0.865265   -0.900453   -0.937082    -0.97521     -1.0149    -1.05621    -1.09921    -1.14398    -1.19057    -1.23907    -1.28955    -1.34209    -1.39678    -1.45371    -1.51296    -1.57464    -1.63883    -1.70565    -1.77519    -1.84758    -1.92292    -2.00134

Independent variable at boundary nodes (tBC) = 0 4

State vector at boundary nodes (xBC) = 
         0         -2
-0.0732871   -2.00134

====================================================================================================

========================
Boundary Value Problem 2
========================

Boundary nodes (tBC) = 
  0 0.5   1

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
   0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09  0.1 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19  0.2 0.21 0.22 0.23 0.24 0.25 0.26 0.27 0.28 0.29  0.3 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39  0.4 0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49  0.5 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59  0.6 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69  0.7 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79  0.8 0.81 0.82 0.83 0.84 0.85 0.86 0.87 0.88 0.89  0.9 0.91 0.92 0.93 0.94 0.95 0.96 0.97 0.98 0.99    1

State vector (x) = 
  -0.0121071   -0.0120974   -0.0120689   -0.0120226   -0.0119594     -0.01188   -0.0117854   -0.0116762   -0.0115532   -0.0114171   -0.0112685   -0.0111081   -0.0109365   -0.0107541   -0.0105617   -0.0103597   -0.0101486  -0.00992884  -0.00970091  -0.00946523  -0.00922221  -0.00897224  -0.00871571  -0.00845296  -0.00818434  -0.00791018   -0.0076308   -0.0073465  -0.00705755  -0.00676426  -0.00646687  -0.00616565  -0.00586084  -0.00555268   -0.0052414  -0.00492723  -0.00461037  -0.00429104  -0.00396944  -0.00364576  -0.00332019  -0.00299293  -0.00266414  -0.00233402  -0.00200273  -0.00167045  -0.00133734  -0.00100358 -0.000669325 -0.000334744 -8.77731e-19  0.000334744  0.000669325   0.00100358   0.00133734   0.00167045   0.00200273   0.00233402   0.00266414   0.00299293   0.00332019   0.00364576   0.00396944   0.00429104   0.00461037   0.00492723    0.0052414   0.00555268   0.00586084   0.00616565   0.00646687   0.00676426   0.00705755    0.0073465    0.0076308   0.00791018   0.00818434   0.00845296   0.00871571   0.00897224   0.00922221   0.00946523   0.00970091   0.00992884    0.0101486    0.0103597    0.0105617    0.0107541    0.0109365    0.0111081    0.0112685    0.0114171    0.0115532    0.0116762    0.0117854      0.01188    0.0119594    0.0120226    0.0120689    0.0120974    0.0121071
-4.47821e-26   0.00192404   0.00375287   0.00549107   0.00714297   0.00871271    0.0102042    0.0116212    0.0129673    0.0142457    0.0154598    0.0166124    0.0177067    0.0187451    0.0197304     0.020665    0.0215513    0.0223915    0.0231876    0.0239417    0.0246556    0.0253311      0.02597    0.0265738     0.027144    0.0276821    0.0281894    0.0286671    0.0291165    0.0295387    0.0299347    0.0303056    0.0306522    0.0309754    0.0312761     0.031555    0.0318127      0.03205    0.0322674    0.0324654    0.0326447    0.0328055    0.0329483    0.0330735    0.0331814    0.0332722    0.0333463    0.0334036    0.0334445     0.033469    0.0334772     0.033469    0.0334445    0.0334036    0.0333463    0.0332722    0.0331814    0.0330735    0.0329483    0.0328055    0.0326447    0.0324654    0.0322674      0.03205    0.0318127     0.031555    0.0312761    0.0309754    0.0306522    0.0303056    0.0299347    0.0295387    0.0291165    0.0286671    0.0281894    0.0276821     0.027144    0.0265738      0.02597    0.0253311    0.0246556    0.0239417    0.0231876    0.0223915    0.0215513     0.020665    0.0197304    0.0187451    0.0177067    0.0166124    0.0154598    0.0142457    0.0129673    0.0116212    0.0102042   0.00871271   0.00714297   0.00549107   0.00375287   0.00192404  5.78964e-17
    0.197323     0.187565     0.178277     0.169434     0.161015     0.152999     0.145365     0.138095      0.13117     0.124573     0.118287     0.112298     0.106589     0.101146    0.0959568    0.0910073    0.0862853    0.0817791    0.0774774    0.0733694    0.0694448     0.065694    0.0621074    0.0586761    0.0553915    0.0522454    0.0492299    0.0463376    0.0435611    0.0408936    0.0383283    0.0358589    0.0334791     0.031183    0.0289649    0.0268192    0.0247406    0.0227239     0.020764     0.018856    0.0169951    0.0151768    0.0133964    0.0116495    0.0099317   0.00823876   0.00656642    0.0049105   0.00326686   0.00163139  7.19398e-12  -0.00163139  -0.00326686   -0.0049105  -0.00656642  -0.00823876   -0.0099317   -0.0116495   -0.0133964   -0.0151768   -0.0169951    -0.018856    -0.020764   -0.0227239   -0.0247406   -0.0268192   -0.0289649    -0.031183   -0.0334791   -0.0358589   -0.0383283   -0.0408936   -0.0435611   -0.0463376   -0.0492299   -0.0522454   -0.0553915   -0.0586761   -0.0621074    -0.065694   -0.0694448   -0.0733694   -0.0774774   -0.0817791   -0.0862853   -0.0910073   -0.0959568    -0.101146    -0.106589    -0.112298    -0.118287    -0.124573     -0.13117    -0.138095    -0.145365    -0.152999    -0.161015    -0.169434    -0.178277    -0.187565    -0.197323

Independent variable at boundary nodes (tBC) =   0 0.5   1

State vector at boundary nodes (xBC) = 
  -0.0121071 -8.77731e-19    0.0121071
-4.47821e-26    0.0334772  5.78964e-17
    0.197323  7.19398e-12    -0.197323

====================================================================================================

========================
Boundary Value Problem 3
========================

Boundary nodes (tBC) = 
       0 0.523599   1.0472   1.5708   2.0944  3.14159

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
        0 0.0261799 0.0523599 0.0785398   0.10472    0.1309   0.15708   0.18326   0.20944  0.235619  0.261799  0.287979  0.314159  0.340339  0.366519  0.392699  0.418879  0.445059  0.471239  0.497419  0.523599  0.523599  0.549779  0.575959  0.602139  0.628319  0.654498  0.680678  0.706858  0.733038  0.759218  0.785398  0.811578  0.837758  0.863938  0.890118  0.916298  0.942478  0.968658  0.994838   1.02102    1.0472    1.0472   1.07338   1.09956   1.12574   1.15192    1.1781   1.20428   1.23046   1.25664   1.28282     1.309   1.33518   1.36136   1.38754   1.41372    1.4399   1.46608   1.49226   1.51844   1.54462    1.5708    1.5708   1.59698   1.62316   1.64934   1.67552    1.7017   1.72788   1.75406   1.78024   1.80642    1.8326   1.85878   1.88496   1.91114   1.93732    1.9635   1.98968   2.01586   2.04204   2.06822    2.0944    2.0944   2.12058   2.14675   2.17293   2.19911   2.22529   2.25147   2.27765   2.30383   2.33001   2.35619   2.38237   2.40855   2.43473   2.46091   2.48709   2.51327   2.53945   2.56563   2.59181   2.61799   2.64417   2.67035   2.69653   2.72271   2.74889   2.77507   2.80125   2.82743   2.85361   2.87979   2.90597   2.93215   2.95833   2.98451   3.01069   3.03687   3.06305   3.08923   3.11541   3.14159

State vector (x) = 
 1.40454e-28    0.0261769     0.052336    0.0784591     0.104528     0.130526     0.156434     0.182236     0.207912     0.233445     0.258819     0.284015     0.309017     0.333807     0.358368     0.382683     0.406737     0.430511      0.45399     0.477159          0.5         -0.5    -0.477159     -0.45399    -0.430511    -0.406737    -0.382683    -0.358368    -0.333807    -0.309017    -0.284015    -0.258819    -0.233445    -0.207912    -0.182236    -0.156434    -0.130526    -0.104528   -0.0784591    -0.052336   -0.0261769 -1.42748e-12 -1.42748e-12    0.0523539     0.104672     0.156918     0.209057     0.261052     0.312869     0.364471     0.415823     0.466891     0.517638     0.568031     0.618034     0.667614     0.716736     0.765367     0.813473     0.861022     0.907981     0.954318            1 -4.25982e-12    0.0261769     0.052336    0.0784591     0.104528     0.130526     0.156434     0.182236     0.207912     0.233445     0.258819     0.284015     0.309017     0.333807     0.358368     0.382683     0.406737     0.430511      0.45399     0.477159          0.5          0.5     0.522499     0.544639     0.566406     0.587785     0.608761      0.62932     0.649448     0.669131     0.688355     0.707107     0.725374     0.743145     0.760406     0.777146     0.793353     0.809017     0.824126     0.838671      0.85264     0.866025     0.878817     0.891007     0.902585     0.913545      0.92388      0.93358     0.942641     0.951057      0.95882     0.965926      0.97237     0.978148     0.983255     0.987688     0.991445     0.994522     0.996917      0.99863     0.999657            1
           1     0.999657      0.99863     0.996917     0.994522     0.991445     0.987688     0.983255     0.978148      0.97237     0.965926      0.95882     0.951057     0.942641      0.93358      0.92388     0.913545     0.902585     0.891007     0.878817     0.866025     0.866025     0.878817     0.891007     0.902585     0.913545      0.92388      0.93358     0.942641     0.951057      0.95882     0.965926      0.97237     0.978148     0.983255     0.987688     0.991445     0.994522     0.996917      0.99863     0.999657            1            2      1.99931      1.99726      1.99383      1.98904      1.98289      1.97538      1.96651       1.9563      1.94474      1.93185      1.91764      1.90211      1.88528      1.86716      1.84776      1.82709      1.80517      1.78201      1.75763      1.73205            1     0.999657      0.99863     0.996917     0.994522     0.991445     0.987688     0.983255     0.978148      0.97237     0.965926      0.95882     0.951057     0.942641      0.93358      0.92388     0.913545     0.902585     0.891007     0.878817     0.866025     0.866025      0.85264     0.838671     0.824126     0.809017     0.793353     0.777146     0.760406     0.743145     0.725374     0.707107     0.688355     0.669131     0.649448      0.62932     0.608761     0.587785     0.566406     0.544639     0.522499          0.5     0.477159      0.45399     0.430511     0.406737     0.382683     0.358368     0.333807     0.309017     0.284015     0.258819     0.233445     0.207912     0.182236     0.156434     0.130526     0.104528    0.0784591     0.052336    0.0261769  4.01902e-12

Independent variable at boundary nodes (tBC) =        0 0.523599 0.523599   1.0472   1.0472   1.5708   1.5708   2.0944   2.0944  3.14159

State vector at boundary nodes (xBC) = 
 1.40454e-28          0.5         -0.5 -1.42748e-12 -1.42748e-12            1 -4.25982e-12          0.5          0.5            1
           1     0.866025     0.866025            1            2      1.73205            1     0.866025     0.866025  4.01902e-12

====================================================================================================

========================
Boundary Value Problem 4
========================

Boundary nodes (tBC) = 
       0 0.785398   1.5708  2.35619  3.14159

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
        0 0.0314159 0.0628319 0.0942478  0.125664   0.15708  0.188496  0.219911  0.251327  0.282743  0.314159  0.345575  0.376991  0.408407  0.439823  0.471239  0.502655  0.534071  0.565487  0.596903  0.628319  0.659734   0.69115  0.722566  0.753982  0.785398  0.785398  0.816814   0.84823  0.879646  0.911062  0.942478  0.973894   1.00531   1.03673   1.06814   1.09956   1.13097   1.16239   1.19381   1.22522   1.25664   1.28805   1.31947   1.35088    1.3823   1.41372   1.44513   1.47655   1.50796   1.53938    1.5708    1.5708   1.60221   1.63363   1.66504   1.69646   1.72788   1.75929   1.79071   1.82212   1.85354   1.88496   1.91637   1.94779    1.9792   2.01062   2.04204   2.07345   2.10487   2.13628    2.1677   2.19911   2.23053   2.26195   2.29336   2.32478   2.35619   2.35619   2.38761   2.41903   2.45044   2.48186   2.51327   2.54469   2.57611   2.60752   2.63894   2.67035   2.70177   2.73319    2.7646   2.79602   2.82743   2.85885   2.89027   2.92168    2.9531   2.98451   3.01593   3.04734   3.07876   3.11018   3.14159

State vector (x) = 
 2.16514e-28    0.0314108    0.0627905    0.0941083     0.125333     0.156434     0.187381     0.218143      0.24869     0.278991     0.309017     0.338738     0.368125     0.397148     0.425779      0.45399     0.481754     0.509041     0.535827     0.562083     0.587785     0.612907     0.637424     0.661312     0.684547     0.707107     0.707107     0.728969     0.750111     0.770513     0.790155     0.809017     0.827081     0.844328     0.860742     0.876307     0.891007     0.904827     0.917755     0.929776     0.940881     0.951057     0.960294     0.968583     0.975917     0.982287     0.987688     0.992115     0.995562     0.998027     0.999507            1    -0.257684     -0.26603    -0.274114    -0.281927    -0.289463    -0.296712    -0.303669    -0.310326    -0.316676    -0.322715    -0.328434     -0.33383    -0.338896    -0.343628    -0.348021     -0.35207    -0.355771    -0.359122    -0.362118    -0.364757    -0.367036    -0.368952    -0.370505    -0.371692    -0.372512    -0.372964    -0.372964    -0.373049    -0.372765    -0.372113    -0.371094    -0.369709    -0.367959    -0.365846    -0.363372    -0.360539    -0.357351    -0.353809    -0.349919    -0.345683    -0.341106    -0.336193    -0.330948    -0.325376    -0.319483    -0.313274    -0.306757    -0.299937    -0.292821    -0.285415    -0.277729    -0.269768
           1     0.999507     0.998027     0.995562     0.992115     0.987688     0.982287     0.975917     0.968583     0.960294     0.951057     0.940881     0.929776     0.917755     0.904827     0.891007     0.876307     0.860742     0.844328     0.827081     0.809017     0.790155     0.770513     0.750111     0.728969     0.707107     0.707107     0.684547     0.661312     0.637424     0.612907     0.587785     0.562083     0.535827     0.509041     0.481754      0.45399     0.425779     0.397148     0.368125     0.338738     0.309017     0.278991      0.24869     0.218143     0.187381     0.156434     0.125333    0.0941083    0.0627905    0.0314108 -7.18934e-13    -0.269768     -0.26154    -0.253055     -0.24432    -0.235344    -0.226136    -0.216704    -0.207059    -0.197209    -0.187165    -0.176936    -0.166532    -0.155964    -0.145242    -0.134377    -0.123379    -0.112259    -0.101029   -0.0896985   -0.0782798   -0.0667839   -0.0552221   -0.0436058   -0.0319464   -0.0202555  -0.00854464  -0.00854464   0.00317467    0.0148908    0.0265923    0.0382676     0.049905    0.0614933    0.0730208    0.0844763    0.0958484     0.107126     0.118298     0.129353      0.14028     0.151069     0.161709     0.172189       0.1825      0.19263      0.20257      0.21231     0.221841     0.231153     0.240236     0.249083     0.257684

Independent variable at boundary nodes (tBC) =        0 0.785398 0.785398   1.5708   1.5708  2.35619  2.35619  3.14159

State vector at boundary nodes (xBC) = 
 2.16514e-28     0.707107     0.707107            1    -0.257684    -0.372964    -0.372964    -0.269768
           1     0.707107     0.707107 -7.18934e-13    -0.269768  -0.00854464  -0.00854464     0.257684

====================================================================================================

Program ended...

[shivanandvp@computer nlmpbvp]$
```

## Source of sample user-side program - test_nlmpBVP.cpp:
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
//     Journal of Mathematical Analysis and Applications 69.2 (1979): 359-371.
// [2] Welsh, Wayne, and Takeo Ojika. "Multipoint boundary value problems with discontinuities I. Algorithms and applications."
//     Journal of Computational and Applied Mathematics 6.2 (1980): 133-143.

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
#include <iostream>             // For the cout statements
#include "nlmpbvp.hpp"          // For the boundary value problem solver function declarations
#include <Eigen/MPRealSupport>  // For arbitrary precision computation
using namespace std;            // For cout
using namespace Eigen;          // For matrix and vector data types and operations
using namespace mpfr;           // For arbitrary precision computation
// ===============================

// ================
// Functions dxBydt
// ================
// dxBydt = a function that defines the derivative of a state vector x at t -- (nx1)

/* Boundary Value Problem 1 */
VectorXm<mpreal> dxBydt_BVP1(mpreal t, VectorXm<mpreal> x){ 
    VectorXm<mpreal> dxdt(2);
    dxdt(0) = x(1);
    dxdt(1) = -fabs(x(0));
    return dxdt;
}

/* Boundary Value Problem 2 */
VectorXm<mpreal> dxBydt_BVP2(mpreal t, VectorXm<mpreal> x){ 
    VectorXm<mpreal> dxdt(3);
    dxdt(0) = x(1);
    dxdt(1) = x(2);
    dxdt(2) = 25*x(1) - 1;
    return dxdt;
}

/* Boundary Value Problem 3 */
VectorXm<mpreal> dxBydt_BVP3(mpreal t, VectorXm<mpreal> x){ 
    VectorXm<mpreal> dxdt(2);
    dxdt(0) =  x(1);
    dxdt(1) = -x(0);
    return dxdt;
}

/* Boundary Value Problem 4 */
VectorXm<mpreal> dxBydt_BVP4(mpreal t, VectorXm<mpreal> x){ 
    VectorXm<mpreal> dxdt(2);
    dxdt(0) =  x(1);
    dxdt(1) = -x(0);
    return dxdt;
}
// ================

// ====================
// Functions BCResidues
// ====================
// BCResidues = a function that defines the boundary condition residues at nodal state vectors xBC -- (nx1) 

/* Boundary Value Problem 1 */
VectorXm<mpreal> BCResidues_BVP1(MatrixXm<mpreal> xBC){
    VectorXm<mpreal> residues(2);
    residues(0) = xBC(0,0) - 0;
    residues(1) = xBC(0,1) + 2;
    return residues;
}

/* Boundary Value Problem 2 */
VectorXm<mpreal> BCResidues_BVP2(MatrixXm<mpreal> xBC){
    VectorXm<mpreal> residues(3);
    residues(0) = xBC(1,0) - 0;
    residues(1) = xBC(1,2) - 0;
    residues(2) = xBC(0,1) - 0;
    return residues;
}


// BCResidues     = a function that defines the boundary condition residues at the nodal state vectors... -- (n(m-1)x1)
//                  ... on the left and right side of integration intervals, xBCL and xBCR

/* Boundary Value Problem 3 */
VectorXm<mpreal> BCResidues_BVP3(MatrixXm<mpreal> xBCL, MatrixXm<mpreal> xBCR){
    VectorXm<mpreal> residues(10);
    residues(0) = xBCL(0,0)-0;
    residues(1) = xBCL(1,0)-1;
    residues(2) = xBCR(0,0) - xBCL(0,1) - 1;
    residues(3) = xBCR(1,0) - xBCL(1,1) - 0;
    residues(4) = xBCR(1,1) - xBCL(1,2) + 1;
    residues(5) = xBCR(0,1) - xBCL(0,2) + 0;
    residues(6) = xBCR(0,2) - xBCL(0,3) - 1;
    residues(7) = xBCR(1,2) - xBCL(1,3) - (sqrt(3)-1);
    residues(8) = xBCR(0,3) - xBCL(0,4) - 0;
    residues(9) = xBCR(1,3) - xBCL(1,4) - 0;
    return residues;
}

/* Boundary Value Problem 4 */
VectorXm<mpreal> BCResidues_BVP4(MatrixXm<mpreal> xBCL, MatrixXm<mpreal> xBCR){
    VectorXm<mpreal> residues(8);
    residues(0) = xBCL(0,0)-0;
    residues(1) = xBCL(1,0)-1;
    residues(2) = pow(xBCL(0,1),2)*pow(xBCL(1,1),3) - exp(xBCR(1,1))*pow(xBCR(1,0),2)*pow(xBCL(1,3),2) + pow(xBCR(1,3),2)*pow(xBCR(0,2),2) - 0.1859767072;
    residues(3) = pow(xBCR(0,1),2)*pow(xBCR(0,0),3)*pow(xBCL(0,3),2) + pow(xBCL(1,2),2)*exp(xBCR(1,2)) + pow(xBCL(0,2),2)*pow(xBCR(0,3),2) - 0.1261677772;
    residues(4) = xBCR(0,0) - xBCL(0,1) + 0;
    residues(5) = xBCR(1,0) - xBCL(1,1) + 0;
    residues(6) = xBCR(0,2) - xBCL(0,3) + 0;
    residues(7) = xBCR(1,2) - xBCL(1,3) + 0;
    return residues;
}


// =================
// The main function
// =================
// This is where the program execution begins
int main(
    int argc,   // argc = the number of arguments passed to the program
    char **argv // argv = an array of strings that are passed to  the program 
    ){

    mpreal::set_default_prec(64); // Set the number of digits of precision you want for computations

    cout<<endl;
    cout<<"============================================================================="<<endl;
    cout<<"Test: Non-linear multipoint boundary value problem solver (nlmpBVP, nlmpBVP2)"<<endl;
    cout<<"============================================================================="<<endl;
    cout<<"Copyright shivanandvp (shivanandvp.oss@gmail.com)"<<endl;

    // Variable declarations   

    /* Boundary Value Problem 1 */
    RowVectorXm<mpreal> tBC_BVP1(2);            // t_BC           = row vector of values at which the boundary conditions are specified              -- (1xm)
    VectorXm<mpreal>   oxt1_BVP1(2);            // oxt1           = column vector of the guessed initial state                                       -- (nx1)

    /* Boundary Value Problem 2 */
    RowVectorXm<mpreal> tBC_BVP2(3);            // t_BC           = row vector of values at which the boundary conditions are specified              -- (1xm)
    VectorXm<mpreal>   oxt1_BVP2(3);            // oxt1           = column vector of the guessed initial state                                       -- (nx1)

    /* Boundary Value Problem 3 */
    RowVectorXm<mpreal> tBC_BVP3(6);            // t_BC           = row vector of values at which the boundary conditions are specified              -- (1xm)  
    MatrixXm<mpreal> oxt1_BVP3(2,5);            // oxt1           = matrix of the guessed initial state                                              -- (nx(m-1))

    /* Boundary Value Problem 4 */
    RowVectorXm<mpreal> tBC_BVP4(5);            // t_BC           = row vector of values at which the boundary conditions are specified              -- (1xm)
    MatrixXm<mpreal> oxt1_BVP4(2,4);            // oxt1           = matrix of the guessed initial state                                              -- (nx(m-1))

    BVPSolution<mpreal> bvpSolution_BVP1;       // bvpSolution    = the structure in which the solutions of the boundary value problem will be saved
    IVAMParameters<mpreal> ivamParameters_BVP1; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

    BVPSolution<mpreal> bvpSolution_BVP2;       // bvpSolution    = the structure in which the solutions of the boundary value problem will be saved
    IVAMParameters<mpreal> ivamParameters_BVP2; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

    BVPSolution<mpreal> bvpSolution_BVP3;       // bvpSolution    = the structure in which the solutions of the boundary value problem will be saved
    IVAMParameters<mpreal> ivamParameters_BVP3; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

    BVPSolution<mpreal> bvpSolution_BVP4;       // bvpSolution    = the structure in which the solutions of the boundary value problem will be saved
    IVAMParameters<mpreal> ivamParameters_BVP4; // ivamParameters = parameters for the Initial Value Adjusting Method (IVAM)

    // Variable definitions

    /* Boundary Value Problem 1 */
    tBC_BVP1  << 0.0, 4.0;                     // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
    oxt1_BVP1 <<  1,                           // oxt1 = column vector of the guessed initial state                                        -- (nx1) 
             0;  
    // tBC_BVP1  << 0.0, 4.0;                     // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
    // oxt1_BVP1 << -1,                           // oxt1 = column vector of the guessed initial state                                        -- (nx1) 
    //          0;

    /* Boundary Value Problem 2 */
    tBC_BVP2  << 0.0, 0.5, 1.0;                // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
    oxt1_BVP2 << 1,                            // oxt1 = column vector of the guessed initial state                                        -- (nx1)
            1,
            1;

    /* Boundary Value Problem 3 */
    // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
    tBC_BVP3  << 0.0, mpfr::const_pi()/6, mpfr::const_pi()/3, mpfr::const_pi()/2, 2*mpfr::const_pi()/3, mpfr::const_pi();
    // oxt1 = a matrix of the guessed initial state on the left side of each integration interval                                        -- (nx(m-1))
    oxt1_BVP3 <<  0.1, 0.1, 0.4, 0.8, 0.9,
            -0.6, 0.1, 0.9, 2.1, 0.8;

    /* Boundary Value Problem 4 */
    // tBC  = the values of the independent variable t at which boundary conditions are defined -- (1xm)
    tBC_BVP4  << 0.0, mpfr::const_pi()/4, mpfr::const_pi()/2, 3*mpfr::const_pi()/4, mpfr::const_pi();
    // oxt1 = a matrix of the guessed initial state on the left side of each integration interval                                        -- (nx(m-1))
    oxt1_BVP4 << -0.1, 0.7, 0.1,-0.7,
             1.1, 0.7,-1.1,-0.7;

    // Assign the parameters for IVAM

    ivamParameters_BVP1.EPSILON    = 1e-10; // EPSILON    = the state perturbation parameter to probe the differential equation system with
    ivamParameters_BVP1.ALPHA      = 1.0;   // ALPHA      = the relaxation factor to scale the adjustment to the initial condition
    ivamParameters_BVP1.SIGMA      = 1e-14; // SIGMA      = the tolerance for error outside which the solver needs to  iterate further. 
    ivamParameters_BVP1.BETA       = 1e-3;  // BETA       = the deflation factor
    ivamParameters_BVP1.printDebug = true;  // printDebug = specify whether debug messages should be output to the console

    ivamParameters_BVP2.EPSILON    = 1e-10; // EPSILON    = the state perturbation parameter to probe the differential equation system with
    ivamParameters_BVP2.ALPHA      = 1.0;   // ALPHA      = the relaxation factor to scale the adjustment to the initial condition
    ivamParameters_BVP2.SIGMA      = 1e-14; // SIGMA      = the tolerance for error outside which the solver needs to  iterate further. 
    ivamParameters_BVP2.BETA       = 1e-3;  // BETA       = the deflation factor
    ivamParameters_BVP2.printDebug = true;  // printDebug = specify whether debug messages should be output to the console

    ivamParameters_BVP3.EPSILON    = 1e-10; // EPSILON    = the state perturbation parameter to probe the differential equation system with
    ivamParameters_BVP3.ALPHA      = 1.0;   // ALPHA      = the relaxation factor to scale the adjustment to the initial condition
    ivamParameters_BVP3.SIGMA      = 1e-14; // SIGMA      = the tolerance for error outside which the solver needs to  iterate further. 
    ivamParameters_BVP3.BETA       = 1e-3;  // BETA       = the deflation factor
    ivamParameters_BVP3.printDebug = true;  // printDebug = specify whether debug messages should be output to the console

    ivamParameters_BVP4.EPSILON    = 1e-10; // EPSILON    = the state perturbation parameter to probe the differential equation system with
    ivamParameters_BVP4.ALPHA      = 1.0;   // ALPHA      = the relaxation factor to scale the adjustment to the initial condition
    ivamParameters_BVP4.SIGMA      = 1e-14; // SIGMA      = the tolerance for error outside which the solver needs to  iterate further. 
    ivamParameters_BVP4.BETA       = 1e-3;  // BETA       = the deflation factor
    ivamParameters_BVP4.printDebug = true;  // printDebug = specify whether debug messages should be output to the console

    /* Boundary Value Problem 1 */
    cout<<endl;
    cout<<"========================"<<endl;
    cout<<"Boundary Value Problem 1"<<endl;
    cout<<"========================"<<endl;
    cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC_BVP1<<endl;
    cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1_BVP1<<endl;
    cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
    bvpSolution_BVP1 = nlmpBVP<mpreal>(2, 2, 101, tBC_BVP1, oxt1_BVP1, dxBydt_BVP1, BCResidues_BVP1, ivamParameters_BVP1);
    cout<<endl;
    cout<<"========"<<endl;
    cout<<"Solution"<<endl;
    cout<<"========"<<endl;
    cout<<"Independent variable (t) = ";
    cout<<endl<<bvpSolution_BVP1.t<<endl;
    cout<<endl;
    cout<<"State vector (x) = "<<endl;
    cout<<bvpSolution_BVP1.x<<endl;
    cout<<endl;    
    cout<<"Independent variable at boundary nodes (tBC) = ";
    cout<<bvpSolution_BVP1.tBC<<endl;
    cout<<endl;
    cout<<"State vector at boundary nodes (xBC) = "<<endl;
    cout<<bvpSolution_BVP1.xBC<<endl;
    cout<<endl; 
    cout<<"===================================================================================================="<<endl; 

     /* Boundary Value Problem 2 */
    cout<<endl;
    cout<<"========================"<<endl;
    cout<<"Boundary Value Problem 2"<<endl;
    cout<<"========================"<<endl;
    cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC_BVP2<<endl;
    cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1_BVP2<<endl;
    cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
    bvpSolution_BVP2 = nlmpBVP<mpreal>(3, 3, 101, tBC_BVP2, oxt1_BVP2, dxBydt_BVP2, BCResidues_BVP2, ivamParameters_BVP2);
    cout<<endl;
    cout<<"========"<<endl;
    cout<<"Solution"<<endl;
    cout<<"========"<<endl;
    cout<<"Independent variable (t) = ";
    cout<<endl<<bvpSolution_BVP2.t<<endl;
    cout<<endl;
    cout<<"State vector (x) = "<<endl;
    cout<<bvpSolution_BVP2.x<<endl;
    cout<<endl;    
    cout<<"Independent variable at boundary nodes (tBC) = ";
    cout<<bvpSolution_BVP2.tBC<<endl;
    cout<<endl;
    cout<<"State vector at boundary nodes (xBC) = "<<endl;
    cout<<bvpSolution_BVP2.xBC<<endl;
    cout<<endl; 
    cout<<"===================================================================================================="<<endl; 

    /* Boundary Value Problem 3 */
    cout<<endl;
    cout<<"========================"<<endl;
    cout<<"Boundary Value Problem 3"<<endl;
    cout<<"========================"<<endl;
    cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC_BVP3<<endl;
    cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1_BVP3<<endl;
    cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
    bvpSolution_BVP3 = nlmpBVP2<mpreal>(2, 6, 12*10+1, tBC_BVP3, oxt1_BVP3, dxBydt_BVP3, BCResidues_BVP3, ivamParameters_BVP3);
    cout<<endl;
    cout<<"========"<<endl;
    cout<<"Solution"<<endl;
    cout<<"========"<<endl;
    cout<<"Independent variable (t) = ";
    cout<<endl<<bvpSolution_BVP3.t<<endl;
    cout<<endl;
    cout<<"State vector (x) = "<<endl;
    cout<<bvpSolution_BVP3.x<<endl;
    cout<<endl;    
    cout<<"Independent variable at boundary nodes (tBC) = ";
    cout<<bvpSolution_BVP3.tBC<<endl;
    cout<<endl;
    cout<<"State vector at boundary nodes (xBC) = "<<endl;
    cout<<bvpSolution_BVP3.xBC<<endl;
    cout<<endl; 
    cout<<"===================================================================================================="<<endl; 

    /* Boundary Value Problem 4 */
    cout<<endl;
    cout<<"========================"<<endl;
    cout<<"Boundary Value Problem 4"<<endl;
    cout<<"========================"<<endl;
    cout<<endl<<"Boundary nodes (tBC) = "<<endl<<tBC_BVP4<<endl;
    cout<<endl<<"Starting state vector (oxt1) = "<<endl<<oxt1_BVP4<<endl;
    cout<<endl<<"Initiating the BVP solver..."<<endl<<endl;
    bvpSolution_BVP4 = nlmpBVP2<mpreal>(2, 5, 10*10+1, tBC_BVP4, oxt1_BVP4, dxBydt_BVP4, BCResidues_BVP4, ivamParameters_BVP4);
    cout<<endl;
    cout<<"========"<<endl;
    cout<<"Solution"<<endl;
    cout<<"========"<<endl;
    cout<<"Independent variable (t) = ";
    cout<<endl<<bvpSolution_BVP4.t<<endl;
    cout<<endl;
    cout<<"State vector (x) = "<<endl;
    cout<<bvpSolution_BVP4.x<<endl;
    cout<<endl;    
    cout<<"Independent variable at boundary nodes (tBC) = ";
    cout<<bvpSolution_BVP4.tBC<<endl;
    cout<<endl;
    cout<<"State vector at boundary nodes (xBC) = "<<endl;
    cout<<bvpSolution_BVP4.xBC<<endl;
    cout<<endl; 
    cout<<"===================================================================================================="<<endl; 

    cout<<endl<<"Program ended..."<<endl<<endl;

    return 0;
}
// =================
```
