%%Inputs: 
%% x=radial coordinate of data (n points); 
%% y=data (n points); 
%% yerr=data uncertainty n points); 
%% rho=radial cordinate of the fitted quantity (k points)  

ga_prof=[
0.000000E+00   2.228277E+00   5.294201E-02   7.103847E-07   7.275483E+00       0   0.000000E+00   9.968179E-05  -2.791286E+03
1.000000E-02   2.228211E+00   5.266311E-02   6.875996E-02   6.479358E+00       1   1.783498E-04   0.000000E+00   7.899838E+02
2.000000E-02   2.227989E+00   5.186663E-02   1.295872E-01   5.697103E+00       2   7.135416E-04   4.217307E-05   2.052298E+02
3.000000E-02   2.227579E+00   5.061762E-02   1.825849E-01   4.935337E+00       3   1.606003E-03   1.150174E-05   9.558167E+01
4.000000E-02   2.226944E+00   4.898592E-02   2.278843E-01   4.205195E+00       4   2.856496E-03   2.492045E-05   5.544595E+01
5.000000E-02   2.226053E+00   4.704390E-02   2.656553E-01   3.526383E+00       5   4.465969E-03   7.667830E-06   3.623733E+01
6.000000E-02   2.224869E+00   4.486546E-02   2.961201E-01   2.934741E+00       6   6.435660E-03   4.792394E-06   2.554567E+01
7.000000E-02   2.223361E+00   4.252467E-02   3.195744E-01   2.493114E+00       7   8.766994E-03   8.626309E-06   1.897857E+01
8.000000E-02   2.221493E+00   4.009479E-02   3.364157E-01   2.290003E+00       8   1.146163E-02   2.396197E-06   1.465433E+01
9.000000E-02   2.219231E+00   3.764733E-02   3.471868E-01   2.387093E+00       9   1.452191E-02   3.354676E-06   1.165510E+01
1.000000E-01   2.216541E+00   3.525006E-02   3.533077E-01   2.446342E+00      10   1.795125E-02   4.792394E-06   9.489094E+00
1.100000E-01   2.213391E+00   3.295714E-02   3.570060E-01   2.542388E+00      11   2.175439E-02   4.313154E-06   7.873463E+00
1.200000E-01   2.209744E+00   3.081707E-02   3.587726E-01   2.577000E+00      12   2.593592E-02   3.594295E-06   6.636055E+00
1.300000E-01   2.205567E+00   2.888268E-02   3.599472E-01   2.379627E+00      13   3.049957E-02   1.437718E-06   5.667151E+00
1.400000E-01   2.200827E+00   2.721570E-02   3.604839E-01   2.209545E+00      14   3.544896E-02   3.833915E-06   4.894139E+00
1.500000E-01   2.195489E+00   2.588109E-02   3.600603E-01   2.073487E+00      15   4.078510E-02   2.396197E-07   4.267407E+00
1.600000E-01   2.189520E+00   2.493871E-02   3.584297E-01   1.978483E+00      16   4.650874E-02   1.677338E-06   3.752113E+00
1.700000E-01   2.182884E+00   2.443154E-02   3.554126E-01   1.930600E+00      17   5.262053E-02   2.156577E-06   3.323205E+00
1.800000E-01   2.175549E+00   2.437470E-02   3.508934E-01   1.933346E+00      18   5.912143E-02   1.677338E-06   2.962308E+00
1.900000E-01   2.167479E+00   2.474970E-02   3.448130E-01   1.976641E+00      19   6.601256E-02   1.198098E-06   2.655684E+00
2.000000E-01   2.158642E+00   2.550705E-02   3.371465E-01   2.061057E+00      20   7.329325E-02   2.396197E-07   2.392894E+00
2.100000E-01   2.149003E+00   2.657615E-02   3.279370E-01   2.188253E+00      21   8.095922E-02   1.198098E-06   2.165895E+00
2.200000E-01   2.138527E+00   2.787715E-02   3.173222E-01   2.322224E+00      22   8.900286E-02   4.792394E-07   1.968407E+00
2.300000E-01   2.127182E+00   2.933315E-02   3.057776E-01   2.390358E+00      23   9.741885E-02   1.557528E-06   1.795471E+00
2.400000E-01   2.114932E+00   3.087906E-02   2.936354E-01   2.480775E+00      24   1.062077E-01   7.188590E-07   1.643132E+00
2.500000E-01   2.101744E+00   3.245959E-02   2.810479E-01   2.598957E+00      25   1.153674E-01   1.198098E-06   1.508193E+00
2.600000E-01   2.087584E+00   3.402775E-02   2.682459E-01   2.730991E+00      26   1.248903E-01   5.990492E-07   1.388059E+00
2.700000E-01   2.072418E+00   3.554344E-02   2.555142E-01   2.841045E+00      27   1.347663E-01   1.258003E-06   1.280599E+00
2.800000E-01   2.056211E+00   3.697347E-02   2.432131E-01   2.908940E+00      28   1.449854E-01   4.193344E-07   1.184056E+00
2.900000E-01   2.038931E+00   3.829218E-02   2.317274E-01   2.944317E+00      29   1.555363E-01   1.078289E-06   1.096962E+00
3.000000E-01   2.020542E+00   3.948181E-02   2.215090E-01   2.927735E+00      30   1.664065E-01   2.396197E-07   1.018091E+00
3.100000E-01   2.001010E+00   4.053119E-02   2.128362E-01   2.861648E+00      31   1.775819E-01   4.792394E-07   9.464090E-01
3.200000E-01   1.980302E+00   4.143486E-02   2.057494E-01   2.817047E+00      32   1.890474E-01   3.594295E-07   8.810393E-01
3.300000E-01   1.958410E+00   4.218927E-02   2.001924E-01   2.811041E+00      33   2.007882E-01   4.792394E-07   8.212361E-01
3.400000E-01   1.935378E+00   4.279100E-02   1.961455E-01   2.818230E+00      34   2.127877E-01   5.990492E-07   7.663589E-01
3.500000E-01   1.911259E+00   4.323711E-02   1.936938E-01   2.835631E+00      35   2.250282E-01   2.995246E-07   7.158573E-01
3.600000E-01   1.886104E+00   4.352587E-02   1.930270E-01   2.842714E+00      36   2.374912E-01   6.589541E-07   6.692557E-01
3.700000E-01   1.859965E+00   4.365790E-02   1.943091E-01   2.847165E+00      37   2.501581E-01   0.000000E+00   6.261422E-01
3.800000E-01   1.832893E+00   4.363703E-02   1.975333E-01   2.780566E+00      38   2.630110E-01   5.391443E-07   5.861567E-01
3.900000E-01   1.804940E+00   4.347010E-02   2.025705E-01   2.731307E+00      39   2.760321E-01   2.096672E-07   5.489853E-01
4.000000E-01   1.776159E+00   4.316607E-02   2.090062E-01   2.627116E+00      40   2.892028E-01   2.096672E-07   5.145112E-01
4.100000E-01   1.746600E+00   4.273841E-02   2.160761E-01   2.482734E+00      41   3.025055E-01   3.294771E-07   4.838567E-01
4.200000E-01   1.716316E+00   4.220445E-02   2.232386E-01   2.353708E+00      42   3.159243E-01   8.985738E-08   4.566976E-01
4.300000E-01   1.685357E+00   4.158219E-02   2.301450E-01   2.260100E+00      43   3.294438E-01   0.000000E+00   4.324404E-01
4.400000E-01   1.653777E+00   4.088865E-02   2.365314E-01   2.198422E+00      44   3.430496E-01   5.990492E-08   4.106208E-01
4.500000E-01   1.621626E+00   4.013896E-02   2.422223E-01   2.166244E+00      45   3.567291E-01   2.695722E-07   3.908699E-01
4.600000E-01   1.588956E+00   3.934564E-02   2.471721E-01   2.167192E+00      46   3.704699E-01   8.985738E-08   3.728923E-01
4.700000E-01   1.555819E+00   3.852055E-02   2.512059E-01   2.080457E+00      47   3.842607E-01   3.594295E-07   3.564473E-01
4.800000E-01   1.522267E+00   3.767708E-02   2.540773E-01   2.046895E+00      48   3.980903E-01   2.995246E-08   3.413378E-01
4.900000E-01   1.488352E+00   3.682665E-02   2.557440E-01   2.069095E+00      49   4.119483E-01   2.396197E-07   3.274019E-01
5.000000E-01   1.454124E+00   3.597721E-02   2.562229E-01   2.127410E+00      50   4.258254E-01   2.096672E-07   3.145039E-01
5.100000E-01   1.419637E+00   3.513366E-02   2.556128E-01   2.237922E+00      51   4.397117E-01   4.193344E-07   3.025310E-01
5.200000E-01   1.384941E+00   3.429696E-02   2.541903E-01   2.393925E+00      52   4.535983E-01   2.995246E-08   2.913881E-01
5.300000E-01   1.350088E+00   3.346430E-02   2.523577E-01   2.586127E+00      53   4.674764E-01   0.000000E+00   2.809949E-01
5.400000E-01   1.315130E+00   3.263007E-02   2.505358E-01   2.756984E+00      54   4.813373E-01   1.797148E-07   2.712840E-01
5.500000E-01   1.280119E+00   3.178716E-02   2.491662E-01   2.953989E+00      55   4.951734E-01   0.000000E+00   2.621982E-01
5.600000E-01   1.245106E+00   3.092756E-02   2.487018E-01   3.093541E+00      56   5.089766E-01   1.497623E-07   2.536892E-01
5.700000E-01   1.210144E+00   3.004723E-02   2.491233E-01   3.189378E+00      57   5.227394E-01   0.000000E+00   2.457167E-01
5.800000E-01   1.175284E+00   2.914668E-02   2.503378E-01   3.329372E+00      58   5.364553E-01   1.497623E-07   2.382474E-01
5.900000E-01   1.140577E+00   2.822742E-02   2.527857E-01   3.485182E+00      59   5.501162E-01   1.797148E-07   2.312537E-01
6.000000E-01   1.106075E+00   2.729212E-02   2.561387E-01   3.462731E+00      60   5.637159E-01   1.647385E-07   2.247142E-01
6.100000E-01   1.071831E+00   2.634598E-02   2.596616E-01   3.446483E+00      61   5.772484E-01   2.096672E-07   2.186125E-01
6.200000E-01   1.037895E+00   2.539399E-02   2.629492E-01   3.533431E+00      62   5.907065E-01   1.497623E-07   2.129375E-01
6.300000E-01   1.004320E+00   2.443779E-02   2.661778E-01   3.705511E+00      63   6.040853E-01   1.347861E-07   2.076948E-01
6.400000E-01   9.711568E-01   2.347718E-02   2.694816E-01   3.741220E+00      64   6.173788E-01   1.048336E-07   2.032494E-01
6.500000E-01   9.384577E-01   2.252138E-02   2.715391E-01   3.657415E+00      65   6.305805E-01   5.990492E-08   1.997961E-01
6.600000E-01   9.062865E-01   2.158763E-02   2.720536E-01   3.717055E+00      66   6.436865E-01   1.946910E-07   1.973915E-01
6.700000E-01   8.747430E-01   2.068783E-02   2.709070E-01   3.737562E+00      67   6.566916E-01   1.497623E-08   1.961307E-01
6.800000E-01   8.439325E-01   1.983360E-02   2.673036E-01   3.824970E+00      68   6.695908E-01   1.497623E-08   1.961602E-01
6.900000E-01   8.139607E-01   1.902626E-02   2.613861E-01   4.115566E+00      69   6.823797E-01   2.246435E-07   1.976943E-01
7.000000E-01   7.849331E-01   1.824747E-02   2.541350E-01   4.458278E+00      70   6.950541E-01   1.497623E-08   2.010451E-01
7.100000E-01   7.569551E-01   1.746556E-02   2.465598E-01   4.892707E+00      71   7.076089E-01   1.647385E-07   2.066693E-01
7.200000E-01   7.301323E-01   1.663609E-02   2.410889E-01   5.461175E+00      72   7.200410E-01   1.497623E-08   2.152496E-01
7.300000E-01   7.045705E-01   1.570699E-02   2.406093E-01   5.780494E+00      73   7.323478E-01   2.695722E-07   2.278399E-01
7.400000E-01   6.803752E-01   1.465095E-02   2.456357E-01   6.157089E+00      74   7.445246E-01   1.797148E-07   2.461417E-01
7.500000E-01   6.576519E-01   1.346874E-02   2.570901E-01   6.488207E+00      75   7.565669E-01   3.145009E-07   2.730755E-01
7.600000E-01   6.365061E-01   1.221185E-02   2.697568E-01   6.216053E+00      76   7.684720E-01   2.995246E-07   3.129996E-01
7.700000E-01   6.170434E-01   1.101151E-02   2.743726E-01   6.510650E+00      77   7.802380E-01   2.096672E-07   3.598243E-01
7.800000E-01   5.993695E-01   9.992704E-03   2.707894E-01   6.859772E+00      78   7.918615E-01   2.995246E-07   4.069137E-01
7.900000E-01   5.835897E-01   9.239744E-03   2.603856E-01   6.881195E+00      79   8.033387E-01   2.995246E-07   4.460856E-01
8.000000E-01   5.698096E-01   8.793045E-03   2.438835E-01   6.843188E+00      80   8.146666E-01   1.198098E-07   4.669365E-01
8.100000E-01   5.579516E-01   8.613187E-03   2.245799E-01   6.491597E+00      81   8.258426E-01   2.995246E-07   4.616646E-01
8.200000E-01   5.473785E-01   8.601867E-03   2.059711E-01   6.369561E+00      82   8.368631E-01   2.995246E-07   4.302972E-01
8.300000E-01   5.373469E-01   8.634840E-03   1.896885E-01   6.547271E+00      83   8.477237E-01   2.995246E-07   3.808651E-01
8.400000E-01   5.271128E-01   8.583111E-03   1.784790E-01   7.001414E+00      84   8.584200E-01   2.096672E-07   3.242150E-01
8.500000E-01   5.159329E-01   8.356057E-03   1.737118E-01   7.553298E+00      85   8.689494E-01   0.000000E+00   2.688783E-01
8.600000E-01   5.030635E-01   7.917799E-03   1.798723E-01   8.041594E+00      86   8.793076E-01   7.488115E-08   2.194180E-01
8.700000E-01   4.877611E-01   7.315938E-03   1.934471E-01   7.562038E+00      87   8.894876E-01   2.995246E-08   1.772858E-01
8.800000E-01   4.692819E-01   6.702800E-03   2.035342E-01   7.330407E+00      88   8.994834E-01   8.985738E-08   1.422241E-01
8.900000E-01   4.468843E-01   6.242695E-03   2.226701E-01   7.007850E+00      89   9.092907E-01   3.744058E-08   1.150298E-01
9.000000E-01   4.201633E-01   6.204612E-03   2.469149E-01   7.020913E+00      90   9.189018E-01   8.236927E-08   9.466216E-02
9.100000E-01   3.894286E-01   6.859311E-03   2.806692E-01   6.696140E+00      91   9.283086E-01   0.000000E+00   7.847121E-02
9.200000E-01   3.550819E-01   8.375727E-03   3.183167E-01   4.200924E+00      92   9.374997E-01   5.241681E-08   6.499794E-02
9.300000E-01   3.175243E-01   9.982001E-03   1.435312E-01   3.757161E+01      93   9.464620E-01   1.872029E-08   5.337019E-02
9.400000E-01   2.771573E-01   8.590859E-03   4.499917E-01   7.388362E+01      94   9.551828E-01   7.488115E-09   4.303021E-02
9.500000E-01   2.343822E-01   4.471703E-03   6.829492E-01   4.335434E+01      95   9.636433E-01   4.118463E-08   3.360021E-02
9.600000E-01   1.896003E-01   5.733459E-03   2.826928E-01   3.805599E+01      96   9.717940E-01   2.433637E-08   2.481002E-02
9.700000E-01   1.432127E-01   6.097849E-03   1.791685E-01   2.854195E+01      97   9.796434E-01   9.360144E-09   1.645575E-02
9.800000E-01   9.562128E-02   4.772726E-03   3.656816E-01   1.902797E+01      98   9.869582E-01   9.828152E-09   8.373944E-03
9.900000E-01   4.722719E-02   4.918876E-03   4.990044E-01   9.513983E+00      99   9.936588E-01   8.775135E-11   4.246403E-04
1.000000E+00  -1.568250E-03   8.683330E-03   5.445036E-01   0.000000E+00     100   9.999999E-01   0              0           ];

rho=ga_prof(:,1);
data_fitted=ga_prof(:,2);
err_data_fitted=ga_prof(:,3);

ga_prof=[ 0.0215975	  2.26739	       0.163369       
 0.0228012	  2.28087	 	0.164043	
 0.0405378	  2.29935	 	0.164967	
 0.0564677	  2.34535	 	0.167268	
 0.0568026	  2.20327	 	0.106953	
 0.0590947	  2.11691	       0.0973317	
 0.0616813	  2.36735	 	0.168368	
 0.0972290	  2.29634	 	0.164817	
  0.107262	  2.21765	 	0.160883	
  0.131018	  2.19890	 	0.159945	
  0.137553	  2.28190	 	0.164095	
  0.155268	  2.17481	 	0.158741	
  0.176771	  2.32192	 	0.166096	
  0.180040	  2.16945	 	0.158472	
  0.205361	  2.15177	 	0.157589	
  0.214807	  2.18566	 	0.159283	
  0.218595	  1.95300	       0.0792188	
  0.231321	  2.08993	 	0.154496	
  0.251705	  2.21672	 	0.160836	
  0.257963	  2.03826	 	0.151913	
  0.285342	  1.93795	 	0.146898	
  0.311999	  2.17706	 	0.134573	
  0.312426	  2.11759	 	0.179766	
  0.312426	  2.36704	 	0.185036	
  0.378167	  1.85225	 	0.144188	
  0.378167	  1.98045	 	0.152338	
  0.469413	  1.42646	 	0.105084	
  0.469413	  1.47860	 	0.125568	
  0.518904	  1.46287	 	0.113698	
  0.518904	  1.28526	 	0.119466	
  0.546383	  1.21399	 	0.108712	
  0.546383	  1.23393	 	0.124053	
  0.579937	  1.16485	       0.0559853	
  0.579937	  1.24228	       0.0675301	
  0.612280	  1.10675	       0.0916625	
  0.612280	  1.20086	 	0.106234	
  0.674805	 0.934823	       0.0409455	
  0.674805	 0.817306	       0.0393381	
  0.737804	 0.627412	       0.0239996	
  0.737804	 0.695731	       0.0332812	
  0.776583	 0.609428	       0.0237180	
  0.776583	 0.652358	       0.0278628	
  0.796464	 0.551862	       0.0280288	
  0.796464	 0.611258	       0.0363550	
  0.816743	 0.543947	       0.0247195	
  0.816743	 0.594803	       0.0314794	
  0.833562	 0.529754	       0.0183763	
  0.833562	 0.541885	       0.0230320	
  0.853399	 0.476323	       0.0173586	
  0.853399	 0.493504	       0.0224270	
  0.869635	 0.513454	       0.0300150	
  0.869635	 0.498166	       0.0350508	
  0.877889	 0.537935	       0.0275473	
  0.877889	 0.439462	       0.0284178	
  0.887645	 0.422507	       0.0205964	
  0.887645	 0.458213	       0.0280981	
  0.896122	 0.457818	       0.0234430	
  0.896122	 0.455585	       0.0312022	
  0.906167	 0.417078	       0.0244766	
  0.906167	 0.436751	       0.0323302	
  0.916411	 0.369042	       0.0203437	
  0.916411	 0.453638	       0.0291565	
  0.919494	 0.263193	       0.0631597	
  0.925378	 0.371476	       0.0300463	
  0.925378	 0.384921	       0.0328432	
  0.936104	 0.301336	       0.0246051	
  0.936104	 0.402844	       0.0297516	
  0.945557	 0.253479	       0.0270681	
  0.945557	 0.301013	       0.0309590	
  0.955342	 0.243664	       0.0201958	
  0.955342	 0.179692	       0.0227987	
  0.968975	 0.158673	       0.0173282	
  0.968975	0.0867081	       0.0109558	
  0.978126	 0.142154	       0.0150677	
  0.978126	0.0523724	      0.00732862	
  0.992284	0.0918243	      0.00938258	
  0.992284	0.0470118	      0.00672708];
  
x=ga_prof(:,1);
y=ga_prof(:,2);
yerr=ga_prof(:,3);  	



%Polynomial degree
m=9;
m0=[];%Degree of the polynom to be equal to zero



%% Generate the z variable as a mapping of input x data range into the interval [-1,1]												 
z=((x-min(x))-(max(x)-x))/(max(x)-min(x));									 													 


%%Defining data variables like in manuscript
b=y./yerr; 


%%Building the Vandermonde matrix
A(:,m+1)=ones(length(z),1,class(z));										 
A(:,m)=z;													 
for j = m-1:-1:1												 
   A(:,j)=2*z.*A(:,j+1)-A(:,j+2);  %% recurrence relation  						 
end														 
%Scaling factor coming out of the integral when mapping the independent variable for Chebycheff polynomials	 
%scal=2/(max(x)-min(x));											 
scal=1; 													 


     
%If one or more degree, m0, of polynom have to be zero, then this should enforce a proper treatment
dump=~ismember([1:m+1],m0);
A=A(:,dump);





%%Inverting
[Q,R] = qr(A,0);
a = R\(Q'*b);
a=a(dump); % Getting rid of unwanted degrees





%%Computing chi2 and quantities that might be of interest 

yfit_data=diag(yerr)*A*a;
res=y-yfit_data;         
db=res./yerr;
chisq = sum(db.^2);
deg3dom = max(0,length(y) - (m+2-length(m0)));
normres = norm(res);
C=pinv(A'*A);




%%Computing uncertainties on the coefficients
da=C*A'*b;





%%Fitted profiles on new cordinate and remapping it onto [-1,1]
rho=rho(:);%making sure it's a column
zz=((rho-min(rho))-(max(rho)-rho))/(max(rho)-min(rho));%new cordinate
%zz=2*rho-1;%Faster computation assuming that rho=[0,1];



%%Computing the jacobian
dzzdx=gradient(zz,rho);
dzzdx=2; %This should be equal to 2 in case rho \in [0,1]





%%Constructing the Vandermonde matrix for the new radial cordinate
%%and its gradient
X=[];
W=[];

%%Building the Vandermonde matrix
X(:,m+1)=ones(length(zz),1,class(zz));
W(:,m+1)=zeros(length(zz),1,class(zz));										 
X(:,m)=zz;
W(:,m)=dzzdx;													 
for j=m-1:-1:1												 
   X(:,j)=2*zz.*X(:,j+1)-X(:,j+2);  %% recurrence relation  
   W(:,j)=2*dzzdx.*(X(:,j+1)+zz.*W(:,j+1))-W(:,j+2);						 
end	

yfit=X*a;
yfit_g=gradient(yfit,rho);
yfit_sl=-yfit_g./rho;





%Computing covariance matrices at points (x,x')
cov_gy=X*C*W';									    
var_yy=X*C*X'; 								    
var_gg=W*C*W'; 								    




							    
%%Computing covariance vectors at points (x,x), i.e. taking the diagonal of the matrix										    
sig2p=diag(var_yy);
sig2g=diag(var_gg);
cov=diag(cov_gy);									    




dump1=6*cov.^2./yfit.^4-4*yfit_g.*cov./yfit.^3-24*yfit_g.*sig2p.*cov./yfit_g.^5;			    
dump2=sig2g+yfit_g.^2;									    
dump3=sig2p./yfit.^4+8*sig2p.^2./yfit.^6+(1./yfit+sig2p./yfit.^3+3*sig2p.^2./yfit.^5).^2; 	    
dump4=(-cov./yfit.^2-3*cov.*sig2p./yfit.^4+yfit_g./yfit+yfit_g./yfit.^3.*sig2p+3*yfit_g.*sig2p.^2./yfit.^5).^2;   
dy2=dump1+dump2.*dump3-dump4;							    
dy2simp=(yfit_g./yfit).^2.*(sig2g./yfit_g.^2+sig2p./yfit.^2);				    





%%Monte-Carlo code loop to estimate uncertainties by brute force
N_MC=1000;
n=length(b);
y1_fit=[];
for i=1:N_MC
   b1=b+randn(n,1);
   a1=R\(Q'*b1);
   a1(dump)=a1;
   y1_fit=[y1_fit X*a];
end


keyboard
