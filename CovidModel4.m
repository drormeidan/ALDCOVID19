% COVID model - v.4  5.4.20

%%% This is the code for the simulations accompanying this arXiv paper:
%%% https://arxiv.org/abs/2004.01453
%%% The notation here is consistent with the notation in the paper.

%%% This script calls CM4.m, where the model equations and solution is
%%% implemented

%% CC
clear; fclose all; close all;
clc;

%% Parametrers
% See the paper and CM4.m for further documentation of these variables.
rate=[1.25,1/5,1/4,1/3,1/11,1/11,1/13,1/13]; % Model rates: beta,r_MR,r_SH,r_CV,r_HR,r_HD,r_VR,r_VD

% lognormal distribution parameters
% para1=[log(3.45),log(11.2),log(2.58),log(4.3)];
% para2=[sqrt(2*log(13/11.2)),sqrt(2*log(13/11.2)),sqrt(2*log(5/4.3)),sqrt(2*log(5/4.3))];

% weibull distribution parameters
para1=[1.47,1.47,1.47,1.47];
para2=[4.41,11.04,3.31,5.52];

prob=[0.3,0.78,0.14,0.08,0.85,0.15,0.5,0.5]; % transition probabilities: p_NS,p_M,p_S,p_C,p_HR,p_HD,p_VR,p_VD
vecf=linspace(0,0.6,15); % fraction of defectors. f=0 means full compliance, f=1 is the opposite.
Tshift=7; % length of a shift (in days)
N=10000000; % population size
t0vec=50; % time of intervention (in days)
dt=0.005; % time increment (in days)
betaI=0.75; % beta after intervation
NumberOfFlips=40; % Number of flips between the two cohorts.

%% Real data
% This section contains data doe each simulated country. Source: Johns Hopkins (see SI)

% yData  - Number of deaths per days
% yDataC - Number of confirmed cases per days

% % Italy
% yData =[0 0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	2	3	7	10	12	17	21	29	34	52	79	107	148	197	233	366	463	631	827	827	1266	1441	1809	2158	2503	2978	3405	4032	4825	5476	6077	6820	7503	8215	9134	10023	10779	11591	12428	13155	13915	14681	15362	15887	16523	17127	17669	18279	18849	19468];
% yDataC=[0 0	0	0	0	0	0	0	0	2	2	2	2	2	2	2	3	3	3	3	3	3	3	3	3	3	3	3	3	3	20	62	155	229	322	453	655	888	1128	1694	2036	2502	3089	3858	4636	5883	7375	9172	10149	12462	12462	17660	21157	24747	27980	31506	35713	41035	47021	53578	59138	63927	69176	74386	80589	86498	92472	97689	101739	105792	110574	115242	119827	124632	128948	132547	135586	139422	143626	147577	152271];
% rate(1)=1.5;
% N=60*(10^6);

% % USA
% yData=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	6	7	11	12	14	17	21	22	28	36	40	47	54	63	85	108	118	200	244	307	417	557	706	942	1209	1581	2026	2467	2978	3873	4757	5926	7087	8407	9619	10783	12722	14695	16478	18586	20463];
% yDataC=[1	1	2	2	5	5	5	5	5	7	8	8	11	11	11	11	11	11	11	11	12	12	13	13	13	13	13	13	13	13	15	15	15	51	51	57	58	60	68	74	98	118	149	217	262	402	518	583	959	1281	1663	2179	2727	3499	4632	6421	7783	13747	19273	25600	33276	43847	53740	65778	83836	101657	121465	140909	161831	188172	213372	243762	275586	308853	337072	366667	396223	429052	461437	496535	526396];
% rate(1)=1.4;
% N=328.2*(10^6);

% % Spain
% yData=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	2	3	5	10	17	28	35	54	55	133	195	289	342	533	623	830	1043	1375	1772	2311	2808	3647	4365	5138	5982	6803	7716	8464	9387	10348	11198	11947	12641	13341	14045	14792	15447	16081	16606];
% yDataC=[0	0	0	0	0	0	0	0	0	0	1	1	1	1	1	1	1	1	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	2	6	13	15	32	45	84	120	165	222	259	400	500	673	1073	1695	2277	2277	5232	6391	7798	9942	11748	13910	17963	20410	25374	28768	35136	39885	49515	57786	65719	73235	80110	87956	95923	104118	112065	119199	126168	131646	136675	141942	148220	153222	158273	163027];
% rate(1)=1.9;
% N=46.9*(10^6);

% % Israel
% yData=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	3	5	8	12	12	15	16	20	26	36	40	44	49	57	65	73	86	95	101];
% yDataC=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	1	2	3	4	7	10	10	12	15	20	37	43	61	61	75	79	100	126	155	213	218	250	304	427	529	712	883	1071	1238	2369	2693	3035	3619	4247	4695	5358	6092	6857	7428	7851	8430	8904	9248	9404	9968	10408	10743];
% rate(1)=1.5;
% N=8.9*(10^6);

% % Germany
% yData=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	2	2	3	3	7	9	11	17	24	28	44	67	84	94	123	157	206	267	342	433	533	645	775	920	1107	1275	1444	1584	1810	2016	2349	2607	2767	2736];
% yDataC=[0	0	0	0	0	1	4	4	4	5	8	10	12	12	12	12	13	13	14	14	16	16	16	16	16	16	16	16	16	16	16	16	16	16	17	27	46	48	79	130	159	196	262	482	670	799	1040	1176	1457	1908	2078	3675	4585	5795	7272	9257	12327	15320	19848	22213	24873	29056	32986	37323	43938	50871	57695	62095	66885	71808	77872	84794	91159	96092	100123	103374	107663	113296	118181	122171	124908];
% rate(1)=1.2;
% N=83*(10^6);

% % Norway
% yData=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	3	3	3	3	6	7	7	7	7	10	12	14	14	19	23	25	32	39	44	50	59	62	71	76	89	101	108	113	119];
% yDataC=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	6	15	19	25	32	56	87	108	147	176	205	400	598	702	996	1090	1221	1333	1463	1550	1746	1914	2118	2385	2621	2863	3084	3369	3755	4015	4284	4445	4641	4863	5147	5370	5550	5687	5865	6086	6086	6211	6314	6409];
% rate(1)=0.85;
% N=5.4*(10^6);

% % South Korea
% yData=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	2	2	6	8	10	12	13	13	16	17	28	28	35	35	42	44	50	53	54	60	66	66	72	75	75	81	84	91	94	102	111	111	120	126	131	139	144	152	158	162	165	169	174	177	183	186	192	200	204	208	211];
% yDataC=[1	1	2	2	3	4	4	4	4	11	12	15	15	16	19	23	24	24	25	27	28	28	28	28	28	29	30	31	31	104	204	433	602	833	977	1261	1766	2337	3150	3736	4335	5186	5621	6088	6593	7041	7314	7478	7513	7755	7869	7979	8086	8162	8236	8320	8413	8565	8652	8799	8961	8961	9037	9137	9241	9332	9478	9583	9661	9786	9887	9976	10062	10156	10237	10284	10331	10384	10423	10450	10480];
% rate(1)=0.8;
% N=51.3*(10^6);

% % Colombia
% yData=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	2	3	3	4	6	6	6	10	12	16	17	19	25	32	35	46	50	54	69	80	100];
% yDataC=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	3	9	9	13	22	34	54	65	93	102	128	196	231	277	378	470	491	539	608	702	798	906	1065	1161	1267	1406	1485	1579	1780	2054	2223	2473	2709];
% rate(1)=1.15;
% N=50.8*(10^6);

% % Argentina
% yData=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	1	2	2	2	2	2	2	3	3	4	4	4	6	8	9	13	18	19	23	27	28	36	39	43	44	48	56	63	72	82	83];
% yDataC=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	2	8	12	12	17	19	19	31	34	45	56	68	79	97	128	158	266	301	387	387	502	589	690	745	820	1054	1054	1133	1265	1451	1451	1554	1628	1715	1795	1975	1975];
% rate(1)=1.1;
% N=45.1*(10^6);

% % Belgium
% yData=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	1	2	2	2	2	2	2	3	3	4	4	4	6	8	9	13	18	19	23	27	28	36	39	43	44	48	56	63	72	82	83];
% yDataC=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	2	8	12	12	17	19	19	31	34	45	56	68	79	97	128	158	266	301	387	387	502	589	690	745	820	1054	1054	1133	1265	1451	1451	1554	1628	1715	1795	1975	1975];
% rate(1)=1.3;
% N=11.6*(10^6);

% % Hubei, China
% yData=[17	17	24	40	52	76	125	125	162	204	249	350	414	479	549	618	699	780	871	974	1068	1068	1310	1457	1596	1696	1789	1921	2029	2144	2144	2346	2346	2495	2563	2615	2641	2682	2727	2761	2803	2835	2871	2902	2931	2959	2986	3008	3024	3046	3056	3062	3075	3085	3099	3111	3122	3130	3133	3139	3153	3153	3160	3163	3169	3174	3177	3182	3186	3187	3193	3199	3203	3207	3210	3212	3212	3213	3215	3216	3219];
% yDataC=[444	444	549	761	1058	1423	3554	3554	4903	5806	7153	11177	13522	16678	19665	22112	24953	27100	29631	31728	33366	33366	48206	54406	56249	58182	59989	61682	62031	62442	62662	64084	64084	64287	64786	65187	65596	65914	66337	66907	67103	67217	67332	67466	67592	67666	67707	67743	67760	67773	67781	67786	67790	67794	67798	67799	67800	67800	67800	67800	67800	67800	67801	67801	67801	67801	67801	67801	67801	67801	67802	67802	67802	67803	67803	67803	67803	67803	67803	67803	67803	67803];
% rate(1)=1.1;
% N=58.5*(10^6);

% % NSW
% yData=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	2	2	2	2	2	2	2	2	2	4	5	5	6	6	6	6	7	7	7	7	8	8	8	8	9	10	12	12	16	18	21	21	21	22	23];
% yDataC=[0	0	0	0	3	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	4	6	6	13	22	22	26	28	38	48	55	65	65	92	112	134	171	210	267	307	353	436	669	669	818	1029	1219	1405	1617	1791	2032	2032	2182	2298	2389	2493	2580	2637	2686	2734	2773	2822	2857];
% rate(1)=0.82;
% N=8.1*(10^6);

% % Austria
% yData=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1	1	3	3	4	6	6	8	16	21	28	30	49	58	68	86	108	128	146	158	168	186	204	220	243	273	295	319	337];
% yDataC=[0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	2	2	3	3	9	14	18	21	29	41	55	79	104	131	182	246	302	504	655	860	1018	1332	1646	2013	2388	2814	3582	4474	5283	5588	6909	7657	8271	8788	9618	10180	10711	11129	11524	11781	12051	12297	12639	12942	13244	13555	13806];
% rate(1)=1.1;
% N=8.9*(10^6);
% 

%% Estimate beta - Section A. Panels (a) - (l) in Figure (2)
% the code for fitting beta to the data

% figure(1);
% hold on
% plot(yData/N,'o');
% axes3 = gca(figure(1));
% hold(axes3,'on');
% box(axes3,'on');
% set(axes3,'FontSize',28,'FontWeight','bold','LineWidth',4,...
%     'PlotBoxAspectRatio',[1 1 1]);
% figure(2);
% hold on
% plot(yDataC/N,'o','Color', [166/255,166/255,166/255],'LineWidth',4,'MarkerSize',8);
% axes3 = gca(figure(2));
% hold(axes3,'on');
% box(axes3,'on');
% set(axes3,'FontSize',28,'FontWeight','bold','LineWidth',4,...
%     'PlotBoxAspectRatio',[1 1 1]);
% [yT,tT] = CM4(dt,120,7,0,rate,para1,para2,prob,0,N,1,0); % without intervation
% figure(1)
% plot(tT,yT(:,12)/N);
% figure(2)
% plot(tT,sum(yT(:,6:8)')/N,'Color', [197/255,90/255,17/255],'LineWidth',4,'MarkerSize',8);

%% Estimate beta - Section B. Panel (m) in Figure (2)
% Ebeta=[1.5 1.4 2.3 1.4 1.4 0.95 1 1.3 1.3 1.4 1.3 0.85 1.3];
% histogram(Ebeta,'Normalization','pdf');

%% The impact of Alternating quarantine.- Panels (e)-(g) in Figure (3)
statevec=[0,1]; % see CM4.m for the meaning of the different states
colorMAP=[[0/255,0/255,153/255];[0/255,153/255,153/255]];
for kstate=1:length(statevec)
    state=statevec(kstate);
    for t0idx=1:length(t0vec)
        t0=t0vec(t0idx);
        % [yT,~]=CM4(dt,t0,7,2,0,rate,para1,para2,prob,NumberOfFlips,N,1,betaI,3,0); %Weekend
        [yT,~]=CM4(dt,t0,7,0,0,rate,para1,para2,prob,NumberOfFlips,N,1,betaI,3,0); %Allweek
        Dq=(yT(end,12)/N);
      Hpeak = zeros(length(vecf),1);
      Vpeak = zeros(length(vecf),1);
      Dinf  = zeros(length(vecf),1);
        for kf=1:length(vecf)
            f=vecf(kf);
            [yT,~]=CM4(dt,t0,7,0,f,rate,para1,para2,prob,NumberOfFlips,N,1,betaI,state,0); %Weekend
            [peakH,~]=(max(yT(:,9)/N));
            [peakV,~]=(max(yT(:,10)/N));
            Hpeak(kf)=peakH;
            Vpeak(kf)=peakV; 
            Dinf(kf)=(yT(end,12)/N)-Dq;
        end

        % % (e)-(g)
        figure(2);
        hold on;
        plot(vecf,Dinf,'MarkerSize',12,'Marker','o','LineWidth',4,...
            'LineStyle','none','Color',colorMAP(kstate,:));
        figure(5);
        hold on;
        plot(vecf,Hpeak,'MarkerSize',12,'Marker','o','LineWidth',4,...
            'LineStyle','none','Color',colorMAP(kstate,:));
        plot(vecf,0.003*ones(1,length(vecf)),'LineWidth',12,...
            'LineStyle','--','Color',[0.5,0.5,0.5]);

        figure(6);
        hold on;
        plot(vecf,Vpeak,'MarkerSize',12,'Marker','o','LineWidth',4,...
            'LineStyle','none','Color',colorMAP(kstate,:));
        plot(vecf,0.00075*ones(1,length(vecf)),'LineWidth',12,...
            'LineStyle','--','Color',[0.5,0.5,0.5]);
    end
end
%% Panels (b)-(d) and panel (a)
vecf=[0.1,0.2,0.3];
%statevec=[3,0,1,2,4];  %panel (a)
%colorMAP=[[0/255,51/255,51/255];[0/255,0/255,153/255];[0/255,153/255,153/255];[102/255,0/255,0/255];[197/255,90/255,17/255]]; %panel (a)
statevec=[0,1];
colorMAP=[[0/255,0/255,153/255];[0/255,153/255,153/255];[197/255,90/255,17/255]];

for kf=1:length(vecf)
    f=vecf(kf);
    for kstate=1:length(statevec)
        state=statevec(kstate);
        [yT,tT]=CM4(dt,t0,7,0,f,rate,para1,para2,prob,NumberOfFlips,N,1,betaI,state,0.5);
        figure(6+kf)
        hold on
        plot(tT,sum(yT(:,6:8),2)'/N,'Color',colorMAP(kstate,:),'LineWidth',12)
    end
end
%% Alternating quarantine vs. population-wide quarantine - Figure (6)
RQpvec=[0.5,0.6,0.7,0.75,0.85,1];
colorMAP=[[102/255,0/255,0/255];[255/255,99/255,71/255];[250/255,128/255,114/255];[255/255,140/255,0/255];[255/255,215/255,0/255];[0/255,51/255,51/255];[0/255,0/255,153/255];[197/255,90/255,17/255]];
for kRQp=1:length(RQpvec) 
    RQp=RQpvec(kRQp);
    [yT,tT]=CM4(dt,t0,7,0,f,rate,para1,para2,prob,NumberOfFlips,N,1,betaI,2,RQp);
    figure(20)
    hold on
    plot(tT,sum(yT(:,6:8),2)'/N,'Color',colorMAP(kRQp,:),'LineWidth',4)
end
[yT,tT]=CM4(dt,t0,7,2,f,rate,para1,para2,prob,NumberOfFlips,N,1,betaI,0,0);
plot(tT,sum(yT(:,6:8),2)'/N,'Color',colorMAP(7,:),'LineWidth',4);

%% Edit the figures:
% (e)-(g)
axes2 = gca(figure(2));
hold(axes2,'on');
box(axes2,'on');
set(axes2,'FontSize',43,'FontWeight','bold','LineWidth',10,...
    'PlotBoxAspectRatio',[1 1 1]);
set(gca,'TickLength',[0.1, 0.001]);
xlim([0 0.5])
plots=get(gca, 'Children');
legend(plots(1:2), {'Int. Q','Alt. Q'});

axes1 = gca(figure(5));
hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',43,'FontWeight','bold','LineWidth',10,...
    'PlotBoxAspectRatio',[1 1 1]);
set(gca,'TickLength',[0.1, 0.001]);
xlim([0 0.5])

axes1 = gca(figure(6));
hold(axes1,'on');
box(axes1,'on');
set(axes1,'FontSize',43,'FontWeight','bold','LineWidth',10,...
    'PlotBoxAspectRatio',[1 1 1]);
set(gca,'TickLength',[0.1, 0.001]);
xlim([0 0.5])

% (b)-(d)
axes3 = gca(figure(7));
hold(axes3,'on');
box(axes3,'on');
set(axes3,'FontSize',43,'FontWeight','bold','LineWidth',10,...
    'PlotBoxAspectRatio',[1 1 1]);
set(gca,'TickLength',[0.1, 0.001]);
xlim([50 150])

axes3 = gca(figure(8));
hold(axes3,'on');
box(axes3,'on');
set(axes3,'FontSize',43,'FontWeight','bold','LineWidth',10,...
    'PlotBoxAspectRatio',[1 1 1]);
set(gca,'TickLength',[0.1, 0.001]);
xlim([50 150])

axes3 = gca(figure(9));
hold(axes3,'on');
box(axes3,'on');
set(axes3,'FontSize',43,'FontWeight','bold','LineWidth',10,...
    'PlotBoxAspectRatio',[1 1 1]);
set(gca,'TickLength',[0.1, 0.001]);
xlim([50 150])

axes2 = gca(figure(20));
hold(axes2,'on');
box(axes2,'on');
set(axes2,'FontSize',28,'FontWeight','bold','LineWidth',4,...
    'PlotBoxAspectRatio',[3 4 1]);
plots=get(gca, 'Children');
legend(plots([1,7:-1:2]), {'AQ','50%','60%','70%','75%','80%','100%'});
set(gca,'TickLength',[0.1, 0.001]);
xlim([20 170])