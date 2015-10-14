classdef RadioLocate < handle
    %function [radioloc, cazacobj] = RadioLocate(varargin)
    %dbstop if error
    
    properties
        clusternodes =[];
        NumRefNodes =[];
        cazacobj = [];
        Band   = 11*15e3;
        DF     = 15e3;
        NREf = [];
        freqvec = [];
        Nsymb  = 1023; %Number of time domain OFDM systems
        CP     = 2*161+1;
        GRA=[];
        REarr=[];
        fs     = 7.68e6;
        ts = [];
        Nfft = [];
        fc     = 3.315e6;
        USR    = 1;
        rprime = 23;
        staticdelay = 1e-5;
        Nfft2 = [];
        CP2 = [];
        %radioloc.cazacobj = [];
        ts2  = [];
        %%%%%%%%%%%%%%%%%%%%%%%%
        LI=[];
        UI=[];
        SNR    = 0;
        dopp   = 0;
        frame   = [1];  %0: Chan. Est., 1: Radio Location, 2: Information
        riciankdB= 8;  %Rician K factor in dB
        pdp    = [   -8  , -6  , -4  ,   0  ,   0  ,  -4  ,  -8   ,   -9  ,  -10  ,  -12  ,  -14  ]; %dB scattering response!
        ricipathdel= [0,  0.1,  0.3,  0.5,   0.7,   1.0,   1.3,   15.0,   15.2,   15.7,   17.2,   20.0]*1e-6;  %in seconds
        chanest= [];
        h00t      =  [];
        h00N      =  [];
        h00Nscat  = [];
        h00Nlos =[];
        spatialvec= [5];
        Numref =[];
        pathdel     =  [];
        HQFREQ2  = [];
        HQFREQscat = [];
        HQFREQLos     =[];
        HQFREQLosB    =[];
        HQFREQLosC    =[];
        FREQISI   = [];
        ChanDbg=[];
        ChanDbg2=[];
        h00Nsinc=[];
        SegSize=[];
        NumbSeg=[];
        SimNum=[];
        simode=[];
    end
    
    methods   %constructor method
        
        function obj=RadioLocate(varargin)
            obj.freqvec= [-(obj.NREf-1)/2:(obj.NREf-1)/2];
            obj.NREf     = floor(obj.Band/obj.DF);
            obj.REarr   = zeros(obj.NREf,obj.Nsymb);
            obj.ts           =1/obj.fs;
            obj.Nfft       = obj.fs/obj.DF;
            obj.Nfft2    = obj.Nfft*obj.USR;
            obj.CP2      = obj.CP  *obj.USR;
            obj.ts2        = obj.ts/obj.USR;
            obj.GRA      =zeros(obj.NREf, obj.Nsymb);
            obj.Numref=obj.spatialvec(1);  %Number of nodes (excluding basis) that have a distinct location geography;
            
            for k=1:nargin
                switch k
                    case 1
                        obj.SNR           =varargin{1};
                    case 2
                        obj.riciankdB     =varargin{2};
                    case 3
                        obj.NumRefNodes   =varargin{3};
                    case 4
                        obj.Nsymb         =varargin{4};
                    case 5
                        obj.frame         =varargin{5};
                    case 6
                        obj.pathdel       =varargin{6};
                    case 7
                        TowerIndex        =varargin{7};
                    case 8
                        obj.CP            =varargin{8};
                    case 9
                        obj.fs            =varargin{9};
                        obj.ts            =1/obj.fs;
                        obj.Nfft          =obj.fs/obj.DF;
                        obj.ts2           =obj.ts/obj.USR;
                    case 10
                        obj.Band          =varargin{10};
                    case 11
                        obj.NREf          =varargin{11};
                    case 13
                        obj.SegSize       =varargin{13};
                    case 14
                        obj.NumbSeg       =varargin{14};
                    case 15
                        obj.SimNum        =varargin{15};
                    case 16
                        obj.simode        =varargin{16};
                end
            end
            
            %radioloc.NREf   = floor(radioloc.Band/radioloc.DF);
            switch rem(obj.NREf,2)
                case 0
                    obj.LI=-obj.NREf/2;
                    obj.UI= obj.NREf/2;
                case 1
                    obj.LI=  -(obj.NREf+1)/2;
                    obj.UI=   (obj.NREf-1)/2;
            end
            obj.clusternodes  =obj.NREf-1;
            
            Valcaz.Band    =obj.Band;
            Valcaz.DF      =obj.DF;
            Valcaz.NREf    =obj.NREf;
            Valcaz.P       =obj.rprime;
            Valcaz.Nsymb   =obj.Nsymb;
            Valcaz.CP      =obj.CP;
            Valcaz.fs      =obj.fs;
            Valcaz.Nfft    =obj.fs/obj.DF;
            Valcaz.fc      =obj.fc;
            Valcaz.delay   =obj.staticdelay;
            Valcaz.USR     =obj.USR;
            Valcaz.fs2     =1/obj.ts2;
            Valcaz.ts2     =obj.ts2;
            Valcaz.frametype=obj.frame;
            Valcaz.SNRdB   =obj.SNR;
            Valcaz.SegSize =obj.SegSize;
            Valcaz.NumbSeg =obj.NumbSeg;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            obj.cazacobj=CazacSys(Valcaz);
            obj.GRA=get(obj.cazacobj,'REArr');
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%
            %Prior Approach
            %%%%%%%%%%%
            PathDelays=obj.ricipathdel;
            %PathDelays=obj.pathdel(TowerIndex)+obj.ricipathdel %modified due to
            %channel cancelation on delest
            %%%%%%%%%%%%%%%%%%%%%
            %New Approach
            [obj.h00t] = GenRice(obj.riciankdB,obj.pdp, PathDelays, obj.ts,  obj.dopp);
            h00N0=MultiRes(obj.h00t, PathDelays, obj.ts, 'NonLinear');
            
            %%%%%%%%%%%%%%%%%%%%%
            %Simulation Mode
            %%%%%%%%%%%%%%%%%%%%%
            
            [tmpdat, obj.h00Nsinc, ToDel]=MultiRes([obj.h00t(1), zeros(size(obj.h00t(2:end)))], PathDelays, obj.ts, 'NonLinear');
            obj.h00N = h00N0/sqrt(sum(h00N0.*conj(h00N0)));
            channelgain= sum(obj.h00N.*conj(obj.h00N));
            
            %Scattering
            obj.h00Nscat = [0,obj.h00N(2:end)];
            GLos=sqrt(10^(obj.riciankdB/10)*(sum(obj.h00Nscat.*conj(obj.h00Nscat)))/...
                (obj.h00N(1)*conj(obj.h00N(1))));
            %Line of site
            obj.h00Nlos =  [GLos*obj.h00N(1),zeros(size(obj.h00N(2:end)))];
            obj.h00N(1)=obj.h00Nlos(1); %Added line to adjust h00N with calulated Gain
            %Normalize Channel
            %%%%%%%%%%%%%%%%%%%%%%%%%
            h00N0=obj.h00N;
            obj.h00N = h00N0/sqrt(sum(h00N0.*conj(h00N0)));
%Test       obj.h00N=[-0.898034469086333 + 0.238174733608791i,-0.0504588034522801 - 0.0129364063637937i,-0.0820231323861387 + 0.0427666499077542i,-0.0607485063709230 + 0.0433491348533733i,0.101325088626001 + 0.130446557563630i,-0.154453126179006 - 0.0152935304477881i,-0.0848843473333600 + 0.0220871897571717i,0.0578568245086465 - 0.0359129876297575i,0.0707940720059820 - 0.0626039295039257i,-0.0157904033331037 + 0.0143137050918215i,-0.172762802784991 + 0.0740190169071698i,-0.00174293994144555 + 0.00301678121028476i,0.00215669289709399 - 0.00243685377295228i,-0.00205658501523218 + 0.00192621698985316i,0.00191764441540824 - 0.00155623639252848i,-0.00179116366991295 + 0.00128420729695954i,0.00168195528720537 - 0.00107790683618847i,-0.00158788942252592 + 0.000916787251983428i,0.00150620885006228 - 0.000787711552646124i,-0.00143459161678966 + 0.000682041436841260i,0.00137122691617433 - 0.000593916809187574i,-0.00131471814511800 + 0.000519244228987139i,0.00126397958812880 - 0.000455088199409622i,-0.00121815563960621 + 0.000399293998926354i,0.00117656207024216 - 0.000350247060195370i,-0.00113864398134771 + 0.000306715214294452i,0.00110394557572958 - 0.000267742650503115i,-0.00107208816194861 + 0.000232577007827269i,0.00104275391126921 - 0.000200618209935277i,-0.00101567367711428 + 0.000171381887145972i,0.000990617723595399 - 0.000144472782837960i,-0.000967388568584209 + 0.000119565120217375i,0.000945815386937898 - 9.63879032115686e-05i,-0.000925749581512318 + 7.47137691669296e-05i,0.000907061240212531 - 5.43504345558572e-05i,-0.000889636273860500 + 3.51340584222147e-05i,0.000873374083355680 - 1.69240412153655e-05i,-0.000858185642804467 - 4.01090068477767e-07i,0.000843991912834179 + 1.69469661601856e-05i,-0.000830722518427118 - 3.28060186697075e-05i,0.000818314640485454 + 4.80596239777553e-05i,-0.000806712081460956 - 6.27798551453105e-05i,0.000795864473799643 + 7.70309237218481e-05i,-0.000785726606378105 - 9.08703745002092e-05i,0.000776257849061318 + 0.000104350082169319i,-0.000767421659375141 - 0.000117517088171041i,0.000759185158309754 + 0.000130414307966018i,-0.000751518764665435 - 0.000143081132700645i,0.000744395879267810 + 0.000155553944464438i,-0.000737792611905451 - 0.000167866560599731i,0.000731687545092730 + 0.000180050619603021i,-0.000726061529767381 - 0.000192135918866103i,0.000720897508862992 + 0.000204150712691358i,-0.000716180365390044 - 0.000216121977577996i,0.000711896792223818 + 0.000228075650636442i,-0.000708035181289017 - 0.000240036846076244i,0.000704585530236785 + 0.000252030053992036i,-0.000701539365063084 - 0.000264079325099965i,0.000698889677430970 + 0.000276208444624273i,-0.000696630875728061 - 0.000288441098187718i,0.000694758749139545 + 0.000300801032293996i,-0.000693270444248875 - 0.000313312211799340i,0.000692164453891346 + 0.000325998976645081i,-0.000691440618197527 - 0.000338886200053562i,0.000691100137977983 + 0.000351999450375183i,-0.000691145600811249 - 0.000365365158815202i,0.000691581020433733 + 0.000379010795357474i,-0.000692411890273458 - 0.000392965055352313i,0.000693645252246371 + 0.000407258059440925i,-0.000695289782248569 - 0.000421921569761559i,0.000697355894129638 + 0.000436989225732631i,-0.000699855864355831 - 0.000452496803141580i,0.000702803980059882 + 0.000468482500808287i,-0.000706216713761197 - 0.000484987259754383i,0.000710112928749217 + 0.000502055120618482i,-0.000714514119974574 - 0.000519733626051834i,0.000719444696338572 + 0.000538074276043473i,-0.000724932311559118 - 0.000557133045614934i,0.000731008252377980 + 0.000576970976162669i,-0.000737707894855451 - 0.000597654853994351i,0.000745071241985047 + 0.000619257992419872i,-0.000753143558986494 - 0.000641861137270306i,0.000761976126620245 + 0.000665553520114528i,-0.000771627137948048 - 0.000690434088991270i,0.000782162770517824 + 0.000716612953505879i,-0.000793658474448313 - 0.000744213090120805i,0.000806200528008495 + 0.000773372365011796i,-0.000819887926940641 - 0.000804245946815325i,0.000834834693280986 + 0.000837009201118127i,-0.000851172715619568 - 0.000871861184257106i,0.000869055268269366 + 0.000909028888181910i,-0.000888661405521547 - 0.000948772434027028i,0.000910201494729169 + 0.000991391474309211i,-0.000933924246828922 - 0.00103723314913876i,0.000960125737929425 + 0.00108670206061556i,-0.000989161110544439 - 0.00114027289690810i,0.00102145992899238 + 0.00119850657672587i,-0.00105754659017071 - 0.00126207113243206i,0.00109806783966603 + 0.00133176906401511i,-0.00114383045021480 - 0.00140857367129112i,0.00119585371895241 + 0.00149367806655055i,-0.00125544404567951 - 0.00158856245706012i,0.00132430322269028 + 0.00169508835049627i,-0.00140468962933987 - 0.00181563346718631i,0.00149966510617677 + 0.00195329004726449i,-0.00161348572405959 - 0.00211216534301871i,0.00175224471036252 + 0.00229785364083623i,-0.00192498009891128 - 0.00251821053616628i,0.00214569264122868 + 0.00278469225417753i,-0.00243728591819701 - 0.00311483202429201i,0.00283997269720665 + 0.00353723093483569i,-0.00343145552321095 - 0.00410285264523713i,0.00438413685840262 + 0.00491516667446166i,-0.00617480465606646 - 0.00623411678336626i,0.0108186304140732 + 0.00905271822894907i,-0.0593733956008914 - 0.0317180578100735i,-0.00983306520189889 + 0.00491020150947775i,0.0199424167650221 + 0.0163872898988712i,-0.00749204107803259 + 0.00197519523939653i,0.00581492895829790 - 0.00918704312495299i,-0.00748490434161264 + 0.0323252135499569i,-0.00346827896328312 + 0.0494783419972939i,-0.000253404748851419 - 0.0157828808779227i,0.000665992221433547 + 0.00967589243549627i,-0.000734312805957502 - 0.00702844806429335i,0.000722572155062833 + 0.00550647889068387i,-0.000693315181936869 - 0.00448488350434689i,0.000669340132094995 + 0.00371584877106540i,-0.000665436708734158 - 0.00306667660974830i,0.000702966934921831 + 0.00242804140985618i,-0.000838292155591556 - 0.00161909791543189i,0.00133093555976416 - 6.25801460322324e-05i,-0.0132603923713611 + 0.0316867646410755i,-0.00132023393797041 + 0.00624636432230502i,0.000638697426769401 - 0.00418611243806657i,-0.000461157608954172 + 0.00347069147930637i,0.000400621275685066 - 0.00307996528760803i,-0.000386364322921907 + 0.00282398438761990i,0.000396167341775416 - 0.00264049837354012i,-0.000421295193200955 + 0.00250286294334814i,0.000458113147471217 - 0.00239781403945770i,-0.000505423118971584 + 0.00231822311050907i,0.000563513910460069 - 0.00226028530769993i,-0.000633883942441558 + 0.00222234574217153i,0.000719340003604638 - 0.00220447775911042i,-0.000824444510504033 + 0.00220853728640362i,0.000956472183902177 - 0.00223870130440732i,-0.00112733121054471 + 0.00230277162092416i,0.00135760473352115 - 0.00241506253130764i,-0.00168590748683765 + 0.00260319251954336i,0.00219387025288270 - 0.00292630616940984i,-0.00308880935896855 + 0.00353544675684188i,0.00509643623264519 - 0.00496043570043796i];
%Test            obj.h00N=[-0.198465725267221 + 0.907636748534868i,-0.114193230559517 - 0.00148292220060957i,0.0154055818288963 - 0.116061447070149i,-0.00447147012987856 - 0.0976750868117240i,-0.0247436868347166 + 0.0398605275136317i,0.0562347782751203 - 0.158001271941040i,0.0348447422162557 - 0.0522740107174352i,-0.00782118091409517 - 0.0228280335455761i,0.0282018914757434 - 0.110367647025616i,-0.0115777096138073 + 0.0325555510433069i,0.0739744479153044 + 0.0380646986771298i,-0.00838323179067988 + 0.0128298244870255i,0.00669036454839782 - 0.0102068478698863i,-0.00566567907394997 + 0.00837844979415066i,0.00492819186881411 - 0.00707425702554738i,-0.00435917879606408 + 0.00610490077948414i,0.00390194551862036 - 0.00535911121223265i,-0.00352414746499070 + 0.00476912382704570i,0.00320541459591230 - 0.00429168803575373i,-0.00293206603735978 + 0.00389802009718515i,0.00269448050282152 - 0.00356827126054630i,-0.00248565016689665 + 0.00328833608123752i,0.00230032650353433 - 0.00304792993683608i,-0.00213448780243953 + 0.00283938641784468i,0.00198499279879150 - 0.00265687985596684i,-0.00184934743659284 + 0.00249590804496841i,0.00172554322083471 - 0.00235293928408065i,-0.00161194240940250 + 0.00222516611486738i,0.00150719473305430 - 0.00211033005837002i,-0.00141017586177424 + 0.00200659464853031i,0.00131994119355941 - 0.00191245197143463i,-0.00123569064397712 + 0.00182665286610955i,0.00115674146672451 - 0.00174815410628408i,-0.00108250702551055 + 0.00167607794886023i,0.00101248003602790 - 0.00160968081037622i,-0.000946219206778262 + 0.00154832876449141i,0.000883338493202264 - 0.00149147819480882i,-0.000823498381723695 + 0.00143866038525522i,0.000766398765386323 - 0.00138946914734159i,-0.000711773078214083 + 0.00134355081097365i,0.000659383432989328 - 0.00130059607038646i,-0.000609016564821517 + 0.00126033329769750i,0.000560480426205330 - 0.00122252302615090i,-0.000513601312120093 + 0.00118695337211725i,0.000468221418847851 - 0.00115343621545578i,-0.000424196759566504 + 0.00112180399631796i,0.000381395374825780 - 0.00109190701596461i,-0.000339695787792998 + 0.00106361115196572i,0.000298985663430019 - 0.00103679591589566i,-0.000259160638108009 + 0.00101135279552673i,0.000220123292019032 - 0.000987183834488983i,-0.000181782241423413 + 0.000964200411051199i,0.000144051331533132 - 0.000942322184613363i,-0.000106848913863955 + 0.000921476184069470i,7.00971943357348e-05 - 0.000901596016680149i,-3.37216403810491e-05 + 0.000882621179742392i,-2.34956308060498e-06 - 0.000864496460303601i,3.81860176601109e-05 + 0.000847171410598995i,-7.38555833378680e-05 - 0.000830599888894200i,0.000109424825147097 + 0.000814739657067964i,-0.000144959441102114 - 0.000799552027648377i,0.000180524677311548 + 0.000785001554172739i,-0.000216185735903143 - 0.000771055759710930i,0.000252008181205640 + 0.000757684899218234i,-0.000288058349577900 - 0.000744861752095302i,0.000324403768346307 + 0.000732561441941070i,-0.000361113589506593 - 0.000720761281028884i,0.000398259044175746 + 0.000709440637514030i,-0.000435913924254044 - 0.000698580823819839i,0.000474155098391124 + 0.000688165005062572i,-0.000513063070171608 - 0.000678178126763120i,0.000552722587470689 + 0.000668606861489909i,-0.000593223313221950 - 0.000659439574472591i,0.000634660569439060 + 0.000650666308651452i,-0.000677136168304992 - 0.000642278790093264i,0.000720759346573264 + 0.000634270455226965i,-0.000765647822524208 - 0.000626636501958996i,0.000811928998421800 + 0.000619373967450528i,-0.000859741336008923 - 0.000612481836204995i,0.000909235938289692 + 0.000605961183179720i,-0.000960578377985644 - 0.000599815357963626i,0.00101395082201789 + 0.000594050217723689i,-0.00106955451268992 - 0.000588674418744329i,0.00112761268063392 + 0.000583699779094223i,-0.00118837398298212 - 0.000579141728468778i,0.00125211658392379 + 0.000575019865833003i,-0.00131915302555419 - 0.000571358651524421i,0.00138983607714007 + 0.000568188268468664i,-0.00146456580397953 - 0.000565545697878732i,0.00154379816766548 + 0.000563476069285479i,-0.00162805556452711 - 0.000562034364509234i,0.00171793983805570 + 0.000561287582423435i,-0.00181414847844828 - 0.000561317509345185i,0.00191749496911470 + 0.000562224293444766i,-0.00202893458782326 - 0.000564131098061048i,0.00214959746769314 + 0.000567190219534761i,-0.00228083144618853 - 0.000571591217795480i,0.00242425829900841 + 0.000577571850619508i,-0.00258184856587185 - 0.000585432970890653i,0.00275602265210163 + 0.000595559116180250i,-0.00294978978961189 - 0.000608447420279861i,0.00316694274159694 + 0.000624748931451361i,-0.00341233661526894 - 0.000645328834896341i,0.00369229816107134 + 0.000671356193673455i,-0.00401524408376398 - 0.000704441077473156i,0.00439264678849535 + 0.000746850215628002i,-0.00484060333907095 - 0.000801857626183454i,0.00538250736652653 + 0.000874337445047480i,-0.00605386815231539 - 0.000971814222139219i,0.00691164659920689 + 0.00110643304249456i,-0.00805407633318439 - 0.00129893026535799i,0.00966819725127196 + 0.00158742675701418i,-0.0121651235790209 - 0.00204963512656509i,0.0166805793599876 + 0.00287161665854797i,-0.0280663658417283 - 0.00466240308977987i,0.147548432981329 + 0.0151015143068658i,0.0297965212044914 - 0.0141328590385516i,-0.00456149441000803 - 0.0610782334035343i,-0.00261300590877931 + 0.0186790961045950i,0.0130674156578613 - 0.0176907826420974i,-0.0503550917442930 + 0.0348610079216761i,-0.0817911524430185 + 0.0378019540236016i,0.0272526054279878 - 0.00955115141721052i,-0.0173203482471014 + 0.00495708086816139i,0.0130094316123784 - 0.00318690247411195i,-0.0105552363205802 + 0.00228107415624759i,0.00895841719639006 - 0.00173791065298478i,-0.00783903368989134 + 0.00137153784899788i,0.00702533664563703 - 0.00109371486794537i,-0.00644111932343452 + 0.000846829531372837i,0.00608960978194102 - 0.000559604040811567i,-0.00621480316091920 - 7.65998953709814e-06i,0.0187337124780871 + 0.0105199624035261i,-0.00273505663270025 + 0.00208924746415155i,0.00322551633379999 - 0.00140941277951126i,-0.00320575234145032 + 0.00117842434530000i,0.00308677360013764 - 0.00105669831733310i,-0.00294022707312032 + 0.000980777644899300i,0.00278655937856437 - 0.000929795017123035i,-0.00263248615721562 + 0.000894828016397575i,0.00247984454908901 - 0.000871492917074465i,-0.00232834335167017 + 0.000857531385155090i,0.00217651333090104 - 0.000851883412092855i,-0.00202197432851287 + 0.000854319918679286i,0.00186134188123158 - 0.000865347402639285i,-0.00168981901828787 + 0.000886307896879405i,0.00150033990180965 - 0.000919708872330437i,-0.00128186901835564 + 0.000969948114216153i,0.00101584324583728 - 0.00104488011678748i,-0.000667955537608106 + 0.00115947357029567i,0.000166237497183778 - 0.00134560343398372i,0.000671436728849759 + 0.00168447525789175i,-0.00248169961373706 - 0.00246067925029420i];
            channelgain= sum(obj.h00N.*conj(obj.h00N));
            %%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(obj.simode,'AWGN')==true;
                h00t_original_length=length(h00N0);
                obj.h00N=zeros(1,h00t_original_length);
                obj.h00Nlos=zeros(1,h00t_original_length);
                obj.h00Nscat=zeros(1,h00t_original_length);
                obj.h00N(1)=1;
                obj.h00Nlos(1)=1;
            end
        end
    end
end
