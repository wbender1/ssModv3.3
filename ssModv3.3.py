#a computer program that tries to figure out the most efficacious RSV vaccine for newborns; simulation prime and boost with different antigens; runs on python 3.11
import myfunc_ssModv3_3 as func
import random
import config as params

#####################
#Table of Contents. #
#####################
'''
1. Define Antigens
2. Create Antigen Strains
3. Create Virus Strains
4. Add Antigens & Viruses to List
5. Prime System
6. Define Equations (rates, counts, time, etc)
7. Output File Naming
8. Run Simulation
'''
###################
#Define Antigens. #
###################
#each strain has two antigens AG1 (G epitopes) AG2 (preF epitopes), ep1-4 = G, ep5-10 = preF
#strain 1 8/60 (subtype B)
AG1_ep1_strain1 = func.Epitope("ep1", '31432132442343143244', params.AG1_ep1_immunogenicity, params.AG1_ep1_neutralization)#G_CCD_AAA47413.1_8/60 (subgroup B)_AG1
AG1_ep2_strain1 = func.Epitope("ep2", '33341241214131331144', params.AG1_ep2_immunogenicity, params.AG1_ep2_neutralization)#G_HBD_AAA47413.1_8/60 (subgroup B)_AG1
AG1_ep3_strain1 = func.Epitope("ep3", '31123223313412321213', params.AG1_ep3_immunogenicity, params.AG1_ep3_neutralization)#G_MUC1_AAA47413.1_8/60 (subgroup B)_AG1
AG1_ep4_strain1 = func.Epitope("ep4", '32444334122241441122', params.AG1_ep4_immunogenicity, params.AG1_ep4_neutralization)#G_MUC2_AAA47413.1_8/60 (subgroup B)_AG1
AG2_ep1_strain1 = func.Epitope("ep1", '31323122123234212123', params.AG2_ep1_immunogenicity, params.AG2_ep1_neutralization)#F_AS0_AZQ19676.1_8/60 (subgroup B)_AG2
AG2_ep2_strain1 = func.Epitope("ep2", '34214223242242222311', params.AG2_ep2_immunogenicity, params.AG2_ep2_neutralization)#F_ASI_AZQ19676.1_8/60 (subgroup B)_AG2
AG2_ep3_strain1 = func.Epitope("ep3", '34342113121212341312', params.AG2_ep3_immunogenicity, params.AG2_ep3_neutralization)#F_ASIIl_AZQ19676.1_8/60 (subgroup B)_AG2
AG2_ep4_strain1 = func.Epitope("ep4", '34144441133441244141', params.AG2_ep4_immunogenicity, params.AG2_ep4_neutralization)#F_ASIII_AZQ19676.1_8/60 (subgroup B)_AG2
AG2_ep5_strain1 = func.Epitope("ep5", '31133424413424433343', params.AG2_ep5_immunogenicity, params.AG2_ep5_neutralization)#F_ASIV_AZQ19676.1_8/60 (subgroup B)_AG2
AG2_ep6_strain1 = func.Epitope("ep6", '32233312321133133231', params.AG2_ep6_immunogenicity, params.AG2_ep6_neutralization)#F_ASV_AZQ19676.1_8/60 (subgroup B)_AG2
#strain 2
AG1_ep1_strain2 = func.Epitope("ep1", '31432132444343113414', params.AG1_ep1_immunogenicity, params.AG1_ep1_neutralization)#G_CCD_AAC14901.1_A2_AG1
AG1_ep2_strain2 = func.Epitope("ep2", '32311144114131231434', params.AG1_ep2_immunogenicity, params.AG1_ep2_neutralization)#G_HBD_AAC14901.1_A2_AG1
AG1_ep3_strain2 = func.Epitope("ep3", '34223231442111322222', params.AG1_ep3_immunogenicity, params.AG1_ep3_neutralization)#G_MUC1_AAC14901.1_A2_AG1
AG1_ep4_strain2 = func.Epitope("ep4", '32224334242323444342', params.AG1_ep4_immunogenicity, params.AG1_ep4_neutralization)#G_MUC2_AAC14901.1_A2_AG1
AG2_ep1_strain2 = func.Epitope("ep1", '33321323124234212123', params.AG2_ep1_immunogenicity, params.AG2_ep1_neutralization)#F_AS0_AAC14902.1_A2_AG2
AG2_ep2_strain2 = func.Epitope("ep2", '34214223232242223313', params.AG2_ep2_immunogenicity, params.AG2_ep2_neutralization)#F_ASI_AAC14902.1_A2_AG2
AG2_ep3_strain2 = func.Epitope("ep3", '34342113421212341312', params.AG2_ep3_immunogenicity, params.AG2_ep3_neutralization)#F_ASII_AAC14902.1_A2_AG2
AG2_ep4_strain2 = func.Epitope("ep4", '33144441133441144144', params.AG2_ep4_immunogenicity, params.AG2_ep4_neutralization)#F_ASIII_AAC14902.1_A2_AG2
AG2_ep5_strain2 = func.Epitope("ep5", '31333424413424433343', params.AG2_ep5_immunogenicity, params.AG2_ep5_neutralization)#F_ASIV_AAC14902.1_A2_AG2
AG2_ep6_strain2 = func.Epitope("ep6", '32233312321133134231', params.AG2_ep6_immunogenicity, params.AG2_ep6_neutralization)#F_ASV_AAC14902.1_A2_AG2
#strain 3
AG1_ep1_strain3 = func.Epitope("ep1", '31432132442343143244', params.AG1_ep1_immunogenicity, params.AG1_ep1_neutralization)#G_CCD_AAB82435.1_B1_AG1
AG1_ep2_strain3 = func.Epitope("ep2", '33141241214131331144', params.AG1_ep2_immunogenicity, params.AG1_ep2_neutralization)#G_HBD_AAB82435.1_B1_AG1
AG1_ep3_strain3 = func.Epitope("ep3", '31323223313412321212', params.AG1_ep3_immunogenicity, params.AG1_ep3_neutralization)#G_MUC1_AAB82435.1_B1_AG1
AG1_ep4_strain3 = func.Epitope("ep4", '31344334122243443122', params.AG1_ep4_immunogenicity, params.AG1_ep4_neutralization)#G_MUC2_AAB82435.1_B1_AG1
AG2_ep1_strain3 = func.Epitope("ep1", '31323322123234212123', params.AG2_ep1_immunogenicity, params.AG2_ep1_neutralization)#F_AS0_AAB82436.1_B1_AG2
AG2_ep2_strain3 = func.Epitope("ep2", '34214223242242221311', params.AG2_ep2_immunogenicity, params.AG2_ep2_neutralization)#F_ASI_AAB82436.1_B1_AG2
AG2_ep3_strain3 = func.Epitope("ep3", '34342113121212341312', params.AG2_ep3_immunogenicity, params.AG2_ep3_neutralization)#F_ASII_AAB82436.1_B1_AG2
AG2_ep4_strain3 = func.Epitope("ep4", '34144441133441244141', params.AG2_ep4_immunogenicity, params.AG2_ep4_neutralization)#F_ASIII_AAB82436.1_B1_AG2
AG2_ep5_strain3 = func.Epitope("ep5", '31133424413424433343', params.AG2_ep5_immunogenicity, params.AG2_ep5_neutralization)#F_ASIV_AAB82436.1_B1_AG2
AG2_ep6_strain3 = func.Epitope("ep6", '32233312321133133231', params.AG2_ep6_immunogenicity, params.AG2_ep6_neutralization)#F_ASV_AAB82436.1_B1_AG2
#strain 4
AG1_ep1_strain4 = func.Epitope("ep1", '31432132442343144244', params.AG1_ep1_immunogenicity, params.AG1_ep1_neutralization)#G_CCD_AHA83901.1_BA_AG1
AG1_ep2_strain4 = func.Epitope("ep2", '33341241213131331144', params.AG1_ep2_immunogenicity, params.AG1_ep2_neutralization)#G_HBD_AHA83901.1_BA_AG1
AG1_ep3_strain4 = func.Epitope("ep3", '31323323313432321212', params.AG1_ep3_immunogenicity, params.AG1_ep3_neutralization)#G_MUC1_AHA83901.1_BA_AG1
AG1_ep4_strain4 = func.Epitope("ep4", '33242314124143411122', params.AG1_ep4_immunogenicity, params.AG1_ep4_neutralization)#G_MUC2_AHA83901.1_BA_AG1
AG2_ep1_strain4 = func.Epitope("ep1", '31323322123234212123', params.AG2_ep1_immunogenicity, params.AG2_ep1_neutralization)#F_AS0_AHA83902.1_BA_AG2
AG2_ep2_strain4 = func.Epitope("ep2", '34214223242242222311', params.AG2_ep2_immunogenicity, params.AG2_ep2_neutralization)#F_ASI_AHA83902.1_BA_AG2
AG2_ep3_strain4 = func.Epitope("ep3", '34342113121212341312', params.AG2_ep3_immunogenicity, params.AG2_ep3_neutralization)#F_ASII_AHA83902.1_BA_AG2
AG2_ep4_strain4 = func.Epitope("ep4", '34144443133441144141', params.AG2_ep4_immunogenicity, params.AG2_ep4_neutralization)#F_ASIII_AHA83902.1_BA_AG2
AG2_ep5_strain4 = func.Epitope("ep5", '31133424413424433343', params.AG2_ep5_immunogenicity, params.AG2_ep5_neutralization)#F_ASIV_AHA83902.1_BA_AG2
AG2_ep6_strain4 = func.Epitope("ep6", '32233312221133134231', params.AG2_ep6_immunogenicity, params.AG2_ep6_neutralization)#F_ASV_AHA83902.1_BA_AG2
#strain 5
AG1_ep1_strain5 = func.Epitope("ep1", '31432132444343113414', params.AG1_ep1_immunogenicity, params.AG1_ep1_neutralization)#G_CCD_AGG39465.1_Bernett_AG1
AG1_ep2_strain5 = func.Epitope("ep2", '32311144214431231434', params.AG1_ep2_immunogenicity, params.AG1_ep2_neutralization)#G_HBD_AGG39465.1_Bernett_AG1
AG1_ep3_strain5 = func.Epitope("ep3", '32223431342111322222', params.AG1_ep3_immunogenicity, params.AG1_ep3_neutralization)#G_MUC1_AGG39465.1_Bernett_AG1
AG1_ep4_strain5 = func.Epitope("ep4", '32224332242323441342', params.AG1_ep4_immunogenicity, params.AG1_ep4_neutralization)#G_MUC2_AGG39465.1_Bernett_AG1
AG2_ep1_strain5 = func.Epitope("ep1", '33321323124234212123', params.AG2_ep1_immunogenicity, params.AG2_ep1_neutralization)#F_AS0_AGG39466.1_Bernett_AG2
AG2_ep2_strain5 = func.Epitope("ep2", '34214221232242223311', params.AG2_ep2_immunogenicity, params.AG2_ep2_neutralization)#F_ASI_AGG39466.1_Bernett_AG2
AG2_ep3_strain5 = func.Epitope("ep3", '34342113421212341312', params.AG2_ep3_immunogenicity, params.AG2_ep3_neutralization)#F_ASII_AGG39466.1_Bernett_AG2
AG2_ep4_strain5 = func.Epitope("ep4", '33144441133441144144', params.AG2_ep4_immunogenicity, params.AG2_ep4_neutralization)#F_ASIII_AGG39466.1_Bernett_AG2
AG2_ep5_strain5 = func.Epitope("ep5", '31333424413424433343', params.AG2_ep5_immunogenicity, params.AG2_ep5_neutralization)#F_ASIV_AGG39466.1_Bernett_AG2
AG2_ep6_strain5 = func.Epitope("ep6", '32233312321133134231', params.AG2_ep6_immunogenicity, params.AG2_ep6_neutralization)#F_ASV_AGG39466.1_Bernett_AG2
#strain 6
AG1_ep1_strain6 = func.Epitope("ep1", '31432132444343113414', params.AG1_ep1_immunogenicity, params.AG1_ep1_neutralization)#G_CCD_AHX57525.1_BJ/35882_AG1
AG1_ep2_strain6 = func.Epitope("ep2", '32311144214131243434', params.AG1_ep2_immunogenicity, params.AG1_ep2_neutralization)#G_HBD_AHX57525.1_BJ/35882_AG1
AG1_ep3_strain6 = func.Epitope("ep3", '34221131344111322222', params.AG1_ep3_immunogenicity, params.AG1_ep3_neutralization)#G_MUC1_AHX57525.1_BJ/35882_AG1
AG1_ep4_strain6 = func.Epitope("ep4", '32324434242423421142', params.AG1_ep4_immunogenicity, params.AG1_ep4_neutralization)#G_MUC2_AHX57525.1_BJ/35882_AG1
AG2_ep1_strain6 = func.Epitope("ep1", '33321323124234212123', params.AG2_ep1_immunogenicity, params.AG2_ep1_neutralization)#F_AS0_AHX57526.1_BJ/35882_AG2
AG2_ep2_strain6 = func.Epitope("ep2", '34214222232242223311', params.AG2_ep2_immunogenicity, params.AG2_ep2_neutralization)#F_ASI_AHX57526.1_BJ/35882_AG2
AG2_ep3_strain6 = func.Epitope("ep3", '34342113421212341312', params.AG2_ep3_immunogenicity, params.AG2_ep3_neutralization)#F_ASII_AHX57526.1_BJ/35882_AG2
AG2_ep4_strain6 = func.Epitope("ep4", '33144441133441144144', params.AG2_ep4_immunogenicity, params.AG2_ep4_neutralization)#F_ASIII_AHX57526.1_BJ/35882_AG2
AG2_ep5_strain6 = func.Epitope("ep5", '31333424413424433343', params.AG2_ep5_immunogenicity, params.AG2_ep5_neutralization)#F_ASIV_AHX57526.1_BJ/35882_AG2
AG2_ep6_strain6 = func.Epitope("ep6", '32233312321133134231', params.AG2_ep6_immunogenicity, params.AG2_ep6_neutralization)#F_ASV_AHX57526.1_BJ/35882_AG2
#strain 7
AG1_ep1_strain7 = func.Epitope("ep1", '31432132444343113414', params.AG1_ep1_immunogenicity, params.AG1_ep1_neutralization)#G_CCD_BAO49759.1_BPH-13-025 (ON1)_AG1
AG1_ep2_strain7 = func.Epitope("ep2", '32311144214131231434', params.AG1_ep2_immunogenicity, params.AG1_ep2_neutralization)#G_HBD_BAO49759.1_BPH-13-025 (ON1)_AG1
AG1_ep3_strain7 = func.Epitope("ep3", '33224231344111322222', params.AG1_ep3_immunogenicity, params.AG1_ep3_neutralization)#G_MUC1_BAO49759.1_BPH-13-025 (ON1)_AG1
AG1_ep4_strain7 = func.Epitope("ep4", '32421324213323421243', params.AG1_ep4_immunogenicity, params.AG1_ep4_neutralization)#G_MUC2_BAO49759.1_BPH-13-025 (ON1)_AG1
AG2_ep1_strain7 = func.Epitope("ep1", '33321323124234212123', params.AG2_ep1_immunogenicity, params.AG2_ep1_neutralization)#F_AS0_AWV55646.1_BPH-13-025 (ON1)_AG2
AG2_ep2_strain7 = func.Epitope("ep2", '34214223232242223311', params.AG2_ep2_immunogenicity, params.AG2_ep2_neutralization)#F_ASI_AWV55646.1_BPH-13-025 (ON1)_AG2
AG2_ep3_strain7 = func.Epitope("ep3", '34342113121212341312', params.AG2_ep3_immunogenicity, params.AG2_ep3_neutralization)#F_ASII_AWV55646.1_BPH-13-025 (ON1)_AG2
AG2_ep4_strain7 = func.Epitope("ep4", '33144441133441144144', params.AG2_ep4_immunogenicity, params.AG2_ep4_neutralization)#F_ASIII_AWV55646.1_BPH-13-025 (ON1)_AG2
AG2_ep5_strain7 = func.Epitope("ep5", '31133424413424433343', params.AG2_ep5_immunogenicity, params.AG2_ep5_neutralization)#F_ASIV_AWV55646.1_BPH-13-025 (ON1)_AG2
AG2_ep6_strain7 = func.Epitope("ep6", '32233312321133134231', params.AG2_ep6_immunogenicity, params.AG2_ep6_neutralization)#F_ASV_AWV55646.1_BPH-13-025 (ON1)_AG2
#strain 8
AG1_ep1_strain8 = func.Epitope("ep1", '31432132441343143244', params.AG1_ep1_immunogenicity, params.AG1_ep1_neutralization)#G_CCD_AGG39486.1_CH-18537_AG1
AG1_ep2_strain8 = func.Epitope("ep2", '33341241214131331144', params.AG1_ep2_immunogenicity, params.AG1_ep2_neutralization)#G_HBD_AGG39486.1_CH-18537_AG1
AG1_ep3_strain8 = func.Epitope("ep3", '31323223313412321213', params.AG1_ep3_immunogenicity, params.AG1_ep3_neutralization)#G_MUC1_AGG39486.1_CH-18537_AG1
AG1_ep4_strain8 = func.Epitope("ep4", '31444334122221441122', params.AG1_ep4_immunogenicity, params.AG1_ep4_neutralization)#G_MUC2_AGG39486.1_CH-18537_AG1
AG2_ep1_strain8 = func.Epitope("ep1", '31323322123234212143', params.AG2_ep1_immunogenicity, params.AG2_ep1_neutralization)#F_AS0_AGG39487.1_CH-18537_AG2
AG2_ep2_strain8 = func.Epitope("ep2", '34214223242242222311', params.AG2_ep2_immunogenicity, params.AG2_ep2_neutralization)#F_ASI_AGG39487.1_CH-18537_AG2
AG2_ep3_strain8 = func.Epitope("ep3", '34342113121212341312', params.AG2_ep3_immunogenicity, params.AG2_ep3_neutralization)#F_ASII_AGG39487.1_CH-18537_AG2
AG2_ep4_strain8 = func.Epitope("ep4", '34144441133441244141', params.AG2_ep4_immunogenicity, params.AG2_ep4_neutralization)#F_ASIII_AGG39487.1_CH-18537_AG2
AG2_ep5_strain8 = func.Epitope("ep5", '31133424413424433343', params.AG2_ep5_immunogenicity, params.AG2_ep5_neutralization)#F_ASIV_AGG39487.1_CH-18537_AG2
AG2_ep6_strain8 = func.Epitope("ep6", '32233312321133133231', params.AG2_ep6_immunogenicity, params.AG2_ep6_neutralization)#F_ASV_AGG39487.1_CH-18537_AG2
#strain 9
AG1_ep1_strain9 = func.Epitope("ep1", '31432132442343143244', params.AG1_ep1_immunogenicity, params.AG1_ep1_neutralization)#G_CCD_AGL98336.1_CH09 G optimized #PB1_AG1
AG1_ep2_strain9 = func.Epitope("ep2", '33121241234131331144', params.AG1_ep2_immunogenicity, params.AG1_ep2_neutralization)#G_HBD_AGL98336.1_CH09 G optimized #PB1_AG1
AG1_ep3_strain9 = func.Epitope("ep3", '31323323313412323213', params.AG1_ep3_immunogenicity, params.AG1_ep3_neutralization)#G_MUC1_AGL98336.1_CH09 G optimized #PB1_AG1
AG1_ep4_strain9 = func.Epitope("ep4", '31444334122143443122', params.AG1_ep4_immunogenicity, params.AG1_ep4_neutralization)#G_MUC2_AGL98336.1_CH09 G optimized #PB1_AG1
AG2_ep1_strain9 = func.Epitope("ep1", '31323322123234212123', params.AG2_ep1_immunogenicity, params.AG2_ep1_neutralization)#F_AS0_BBD49539.1_PH15 F optimized #PA1_AG2
AG2_ep2_strain9 = func.Epitope("ep2", '34214223342242223311', params.AG2_ep2_immunogenicity, params.AG2_ep2_neutralization)#F_ASI_BBD49539.1_PH15 F optimized #PA1_AG2
AG2_ep3_strain9 = func.Epitope("ep3", '34342113121212241312', params.AG2_ep3_immunogenicity, params.AG2_ep3_neutralization)#F_ASII_BBD49539.1_PH15 F optimized #PA1_AG2
AG2_ep4_strain9 = func.Epitope("ep4", '33144442133441444141', params.AG2_ep4_immunogenicity, params.AG2_ep4_neutralization)#F_ASIII_BBD49539.1_PH15 F optimized #PA1_AG2
AG2_ep5_strain9 = func.Epitope("ep5", '31133424413424432323', params.AG2_ep5_immunogenicity, params.AG2_ep5_neutralization)#F_ASIV_BBD49539.1_PH15 F optimized #PA1_AG2
AG2_ep6_strain9 = func.Epitope("ep6", '32233312321133144231', params.AG2_ep6_immunogenicity, params.AG2_ep6_neutralization)#F_ASV_BBD49539.1_PH15 F optimized #PA1_AG2
#strain 10
AG1_ep1_strain10 = func.Epitope("ep1", '34432132444333113414', params.AG1_ep1_immunogenicity, params.AG1_ep1_neutralization)#G_CCD_QDF83109.1_CH18 G optimized #PB2_AG1
AG1_ep2_strain10 = func.Epitope("ep2", '32311144114131231434', params.AG1_ep2_immunogenicity, params.AG1_ep2_neutralization)#G_HBD_QDF83109.1_CH18 G optimized #PB2_AG1
AG1_ep3_strain10 = func.Epitope("ep3", '33224231144111322222', params.AG1_ep3_immunogenicity, params.AG1_ep3_neutralization)#G_MUC1_QDF83109.1_CH18 G optimized #PB2_AG1
AG1_ep4_strain10 = func.Epitope("ep4", '32421323213323424242', params.AG1_ep4_immunogenicity, params.AG1_ep4_neutralization)#G_MUC2_QDF83109.1_CH18 G optimized #PB2_AG1
AG2_ep1_strain10 = func.Epitope("ep1", '33321323124234212123', params.AG2_ep1_immunogenicity, params.AG2_ep1_neutralization)#F_AS0_BBB35187.1_JA04 F optimized #PA2_AG2
AG2_ep2_strain10 = func.Epitope("ep2", '34214221232242223111', params.AG2_ep2_immunogenicity, params.AG2_ep2_neutralization)#F_ASI_BBB35187.1_JA04 F optimized #PA2_AG2
AG2_ep3_strain10 = func.Epitope("ep3", '34342113421212341312', params.AG2_ep3_immunogenicity, params.AG2_ep3_neutralization)#F_ASII_BBB35187.1_JA04 F optimized #PA2_AG2
AG2_ep4_strain10 = func.Epitope("ep4", '33144441133441144144', params.AG2_ep4_immunogenicity, params.AG2_ep4_neutralization)#F_ASIII_BBB35187.1_JA04 F optimized #PA2_AG2
AG2_ep5_strain10 = func.Epitope("ep5", '31333424413424433343', params.AG2_ep5_immunogenicity, params.AG2_ep5_neutralization)#F_ASIV_BBB35187.1_JA04 F optimized #PA2_AG2
AG2_ep6_strain10 = func.Epitope("ep6", '32233312321133134231', params.AG2_ep6_immunogenicity, params.AG2_ep6_neutralization)#F_ASV_BBB35187.1_JA04 F optimized #PA2_AG2
#strain 11
AG1_ep1_strain11 = func.Epitope("ep1", '31432132444343113414', params.AG1_ep1_immunogenicity, params.AG1_ep1_neutralization)#G_CCD_AGG39393.1_Long_AG1
AG1_ep2_strain11 = func.Epitope("ep2", '32311144414131231434', params.AG1_ep2_immunogenicity, params.AG1_ep2_neutralization)#G_HBD_AGG39393.1_Long_AG1
AG1_ep3_strain11 = func.Epitope("ep3", '34224431343111322222', params.AG1_ep3_immunogenicity, params.AG1_ep3_neutralization)#G_MUC1_AGG39393.1_Long_AG1
AG1_ep4_strain11 = func.Epitope("ep4", '32124234242323441342', params.AG1_ep4_immunogenicity, params.AG1_ep4_neutralization)#G_MUC2_AGG39393.1_Long_AG1
AG2_ep1_strain11 = func.Epitope("ep1", '33321323124234212123', params.AG2_ep1_immunogenicity, params.AG2_ep1_neutralization)#F_AS0_AGG39394.1_Long_AG2
AG2_ep2_strain11 = func.Epitope("ep2", '34214224232242223311', params.AG2_ep2_immunogenicity, params.AG2_ep2_neutralization)#F_ASI_AGG39394.1_Long_AG2
AG2_ep3_strain11 = func.Epitope("ep3", '34342113421212341312', params.AG2_ep3_immunogenicity, params.AG2_ep3_neutralization)#F_ASII_AGG39394.1_Long_AG2
AG2_ep4_strain11 = func.Epitope("ep4", '33144441133441144144', params.AG2_ep4_immunogenicity, params.AG2_ep4_neutralization)#F_ASIII_AGG39394.1_Long_AG2
AG2_ep5_strain11 = func.Epitope("ep5", '31333424413424433343', params.AG2_ep5_immunogenicity, params.AG2_ep5_neutralization)#F_ASIV_AGG39394.1_Long_AG2
AG2_ep6_strain11 = func.Epitope("ep6", '32233312321133134231', params.AG2_ep6_immunogenicity, params.AG2_ep6_neutralization)#F_ASV_AGG39394.1_Long_AG2
#strain 12
AG1_ep1_strain12 = func.Epitope("ep1", '31432132434343113414', params.AG1_ep1_immunogenicity, params.AG1_ep1_neutralization)#G_CCD_AHV82099.1_NA1_AG1
AG1_ep2_strain12 = func.Epitope("ep2", '32311144114131231434', params.AG1_ep2_immunogenicity, params.AG1_ep2_neutralization)#G_HBD_AHV82099.1_NA1_AG1
AG1_ep3_strain12 = func.Epitope("ep3", '33224231344111322222', params.AG1_ep3_immunogenicity, params.AG1_ep3_neutralization)#G_MUC1_AHV82099.1_NA1_AG1
AG1_ep4_strain12 = func.Epitope("ep4", '32424334242323421142', params.AG1_ep4_immunogenicity, params.AG1_ep4_neutralization)#G_MUC2_AHV82099.1_NA1_AG1
AG2_ep1_strain12 = func.Epitope("ep1", '33321323124234212123', params.AG2_ep1_immunogenicity, params.AG2_ep1_neutralization)#F_AS0_AHV82100.1_NA1_AG2
AG2_ep2_strain12 = func.Epitope("ep2", '34214222232242223311', params.AG2_ep2_immunogenicity, params.AG2_ep2_neutralization)#F_ASI_AHV82100.1_NA1_AG2
AG2_ep3_strain12 = func.Epitope("ep3", '34342113421212341312', params.AG2_ep3_immunogenicity, params.AG2_ep3_neutralization)#F_ASII_AHV82100.1_NA1_AG2
AG2_ep4_strain12 = func.Epitope("ep4", '33144441133441144144', params.AG2_ep4_immunogenicity, params.AG2_ep4_neutralization)#F_ASIII_AHV82100.1_NA1_AG2
AG2_ep5_strain12 = func.Epitope("ep5", '31333424413424433343', params.AG2_ep5_immunogenicity, params.AG2_ep5_neutralization)#F_ASIV_AHV82100.1_NA1_AG2
AG2_ep6_strain12 = func.Epitope("ep6", '32233312321133134231', params.AG2_ep6_immunogenicity, params.AG2_ep6_neutralization)#F_ASV_AHV82100.1_NA1_AG2
##########################
#Create Antigen Strains. #
##########################
#add both G and F epitopes to one antigen in list
#strain 1
AG1_strain1 = func.AntigenType("Type0", 0)
AG1_strain1.add_epitope(AG1_ep1_strain1); AG1_strain1.add_epitope(AG1_ep2_strain1); AG1_strain1.add_epitope(AG1_ep3_strain1); AG1_strain1.add_epitope(AG1_ep4_strain1); AG1_strain1.add_epitope(AG2_ep1_strain1); AG1_strain1.add_epitope(AG2_ep2_strain1); AG1_strain1.add_epitope(AG2_ep3_strain1); AG1_strain1.add_epitope(AG2_ep4_strain1); AG1_strain1.add_epitope(AG2_ep5_strain1); AG1_strain1.add_epitope(AG2_ep6_strain1)
#strain 2
AG1_strain2 = func.AntigenType("Type1", 1)
AG1_strain2.add_epitope(AG1_ep1_strain2); AG1_strain2.add_epitope(AG1_ep2_strain2); AG1_strain2.add_epitope(AG1_ep3_strain2); AG1_strain2.add_epitope(AG1_ep4_strain2); AG1_strain2.add_epitope(AG2_ep1_strain2); AG1_strain2.add_epitope(AG2_ep2_strain2); AG1_strain2.add_epitope(AG2_ep3_strain2); AG1_strain2.add_epitope(AG2_ep4_strain2); AG1_strain2.add_epitope(AG2_ep5_strain2); AG1_strain2.add_epitope(AG2_ep6_strain2)
#strain 3
AG1_strain3 = func.AntigenType("Type2", 2)
AG1_strain3.add_epitope(AG1_ep1_strain3); AG1_strain3.add_epitope(AG1_ep2_strain3); AG1_strain3.add_epitope(AG1_ep3_strain3); AG1_strain3.add_epitope(AG1_ep4_strain3); AG1_strain3.add_epitope(AG2_ep1_strain3); AG1_strain3.add_epitope(AG2_ep2_strain3); AG1_strain3.add_epitope(AG2_ep3_strain3); AG1_strain3.add_epitope(AG2_ep4_strain3); AG1_strain3.add_epitope(AG2_ep5_strain3); AG1_strain3.add_epitope(AG2_ep6_strain3)
#strain 4
AG1_strain4 = func.AntigenType("Type3", 3)
AG1_strain4.add_epitope(AG1_ep1_strain4); AG1_strain4.add_epitope(AG1_ep2_strain4); AG1_strain4.add_epitope(AG1_ep3_strain4); AG1_strain4.add_epitope(AG1_ep4_strain4); AG1_strain4.add_epitope(AG2_ep1_strain4); AG1_strain4.add_epitope(AG2_ep2_strain4); AG1_strain4.add_epitope(AG2_ep3_strain4); AG1_strain4.add_epitope(AG2_ep4_strain4); AG1_strain4.add_epitope(AG2_ep5_strain4); AG1_strain4.add_epitope(AG2_ep6_strain4)
#strain 5
AG1_strain5 = func.AntigenType("Type4", 4)
AG1_strain5.add_epitope(AG1_ep1_strain5); AG1_strain5.add_epitope(AG1_ep2_strain5); AG1_strain5.add_epitope(AG1_ep3_strain5); AG1_strain5.add_epitope(AG1_ep4_strain5); AG1_strain5.add_epitope(AG2_ep1_strain5); AG1_strain5.add_epitope(AG2_ep2_strain5); AG1_strain5.add_epitope(AG2_ep3_strain5); AG1_strain5.add_epitope(AG2_ep4_strain5); AG1_strain5.add_epitope(AG2_ep5_strain5); AG1_strain5.add_epitope(AG2_ep6_strain5)
#strain 6
AG1_strain6 = func.AntigenType("Type5", 5)
AG1_strain6.add_epitope(AG1_ep1_strain6); AG1_strain6.add_epitope(AG1_ep2_strain6); AG1_strain6.add_epitope(AG1_ep3_strain6); AG1_strain6.add_epitope(AG1_ep4_strain6); AG1_strain6.add_epitope(AG2_ep1_strain6); AG1_strain6.add_epitope(AG2_ep2_strain6); AG1_strain6.add_epitope(AG2_ep3_strain6); AG1_strain6.add_epitope(AG2_ep4_strain6); AG1_strain6.add_epitope(AG2_ep5_strain6); AG1_strain6.add_epitope(AG2_ep6_strain6)
#strain 7
AG1_strain7 = func.AntigenType("Type6", 6)
AG1_strain7.add_epitope(AG1_ep1_strain7); AG1_strain7.add_epitope(AG1_ep2_strain7); AG1_strain7.add_epitope(AG1_ep3_strain7); AG1_strain7.add_epitope(AG1_ep4_strain7); AG1_strain7.add_epitope(AG2_ep1_strain7); AG1_strain7.add_epitope(AG2_ep2_strain7); AG1_strain7.add_epitope(AG2_ep3_strain7); AG1_strain7.add_epitope(AG2_ep4_strain7); AG1_strain7.add_epitope(AG2_ep5_strain7); AG1_strain7.add_epitope(AG2_ep6_strain7)
#strain 8
AG1_strain8 = func.AntigenType("Type7", 7)
AG1_strain8.add_epitope(AG1_ep1_strain8); AG1_strain8.add_epitope(AG1_ep2_strain8); AG1_strain8.add_epitope(AG1_ep3_strain8); AG1_strain8.add_epitope(AG1_ep4_strain8); AG1_strain8.add_epitope(AG2_ep1_strain8); AG1_strain8.add_epitope(AG2_ep2_strain8); AG1_strain8.add_epitope(AG2_ep3_strain8); AG1_strain8.add_epitope(AG2_ep4_strain8); AG1_strain8.add_epitope(AG2_ep5_strain8); AG1_strain8.add_epitope(AG2_ep6_strain8)
#strain 9
AG1_strain9 = func.AntigenType("Type8", 8)
AG1_strain9.add_epitope(AG1_ep1_strain9); AG1_strain9.add_epitope(AG1_ep2_strain9); AG1_strain9.add_epitope(AG1_ep3_strain9); AG1_strain9.add_epitope(AG1_ep4_strain9); AG1_strain9.add_epitope(AG2_ep1_strain9); AG1_strain9.add_epitope(AG2_ep2_strain9); AG1_strain9.add_epitope(AG2_ep3_strain9); AG1_strain9.add_epitope(AG2_ep4_strain9); AG1_strain9.add_epitope(AG2_ep5_strain9); AG1_strain9.add_epitope(AG2_ep6_strain9)
#strain 10
AG1_strain10 = func.AntigenType("Type9", 9)
AG1_strain10.add_epitope(AG1_ep1_strain10); AG1_strain10.add_epitope(AG1_ep2_strain10); AG1_strain10.add_epitope(AG1_ep3_strain10); AG1_strain10.add_epitope(AG1_ep4_strain10); AG1_strain10.add_epitope(AG2_ep1_strain10); AG1_strain10.add_epitope(AG2_ep2_strain10); AG1_strain10.add_epitope(AG2_ep3_strain10); AG1_strain10.add_epitope(AG2_ep4_strain10); AG1_strain10.add_epitope(AG2_ep5_strain10); AG1_strain10.add_epitope(AG2_ep6_strain10)
#strain 11
AG1_strain11 = func.AntigenType("Type10", 10)
AG1_strain11.add_epitope(AG1_ep1_strain11); AG1_strain11.add_epitope(AG1_ep2_strain11); AG1_strain11.add_epitope(AG1_ep3_strain11); AG1_strain11.add_epitope(AG1_ep4_strain11); AG1_strain11.add_epitope(AG2_ep1_strain11); AG1_strain11.add_epitope(AG2_ep2_strain11); AG1_strain11.add_epitope(AG2_ep3_strain11); AG1_strain11.add_epitope(AG2_ep4_strain11); AG1_strain11.add_epitope(AG2_ep5_strain11); AG1_strain11.add_epitope(AG2_ep6_strain11)
#strain 12
AG1_strain12 = func.AntigenType("Type11", 11)
AG1_strain12.add_epitope(AG1_ep1_strain12); AG1_strain12.add_epitope(AG1_ep2_strain12); AG1_strain12.add_epitope(AG1_ep3_strain12); AG1_strain12.add_epitope(AG1_ep4_strain12); AG1_strain12.add_epitope(AG2_ep1_strain12); AG1_strain12.add_epitope(AG2_ep2_strain12); AG1_strain12.add_epitope(AG2_ep3_strain12); AG1_strain12.add_epitope(AG2_ep4_strain12); AG1_strain12.add_epitope(AG2_ep5_strain12); AG1_strain12.add_epitope(AG2_ep6_strain12)
########################
#Create Virus Strains. #
########################
#add both G and F epitopes to one antigen in list
#virus 1
V1_strain1 = func.AntigenType("Type12", 12)
V1_strain1.add_epitope(AG1_ep1_strain1); V1_strain1.add_epitope(AG1_ep2_strain1); V1_strain1.add_epitope(AG1_ep3_strain1); V1_strain1.add_epitope(AG1_ep4_strain1); V1_strain1.add_epitope(AG2_ep1_strain1); V1_strain1.add_epitope(AG2_ep2_strain1); V1_strain1.add_epitope(AG2_ep3_strain1); V1_strain1.add_epitope(AG2_ep4_strain1); V1_strain1.add_epitope(AG2_ep5_strain1); V1_strain1.add_epitope(AG2_ep6_strain1)
#virus 2
V2_strain2 = func.AntigenType("Type13", 13)
V2_strain2.add_epitope(AG1_ep1_strain2); V2_strain2.add_epitope(AG1_ep2_strain2); V2_strain2.add_epitope(AG1_ep3_strain2); V2_strain2.add_epitope(AG1_ep4_strain2); V2_strain2.add_epitope(AG2_ep1_strain2); V2_strain2.add_epitope(AG2_ep2_strain2); V2_strain2.add_epitope(AG2_ep3_strain2); V2_strain2.add_epitope(AG2_ep4_strain2); V2_strain2.add_epitope(AG2_ep5_strain2); V2_strain2.add_epitope(AG2_ep6_strain2)
#virus 3
V3_strain3 = func.AntigenType("Type14", 14)
V3_strain3.add_epitope(AG1_ep1_strain3); V3_strain3.add_epitope(AG1_ep2_strain3); V3_strain3.add_epitope(AG1_ep3_strain3); V3_strain3.add_epitope(AG1_ep4_strain3); V3_strain3.add_epitope(AG2_ep1_strain3); V3_strain3.add_epitope(AG2_ep2_strain3); V3_strain3.add_epitope(AG2_ep3_strain3); V3_strain3.add_epitope(AG2_ep4_strain3); V3_strain3.add_epitope(AG2_ep5_strain3); V3_strain3.add_epitope(AG2_ep6_strain3)
#virus 4
V4_strain4 = func.AntigenType("Type15", 15)
V4_strain4.add_epitope(AG1_ep1_strain4); V4_strain4.add_epitope(AG1_ep2_strain4); V4_strain4.add_epitope(AG1_ep3_strain4); V4_strain4.add_epitope(AG1_ep4_strain4); V4_strain4.add_epitope(AG2_ep1_strain4); V4_strain4.add_epitope(AG2_ep2_strain4); V4_strain4.add_epitope(AG2_ep3_strain4); V4_strain4.add_epitope(AG2_ep4_strain4); V4_strain4.add_epitope(AG2_ep5_strain4); V4_strain4.add_epitope(AG2_ep6_strain4)
#virus 5
V5_strain5 = func.AntigenType("Type16", 16)
V5_strain5.add_epitope(AG1_ep1_strain5); V5_strain5.add_epitope(AG1_ep2_strain5); V5_strain5.add_epitope(AG1_ep3_strain5); V5_strain5.add_epitope(AG1_ep4_strain5); V5_strain5.add_epitope(AG2_ep1_strain5); V5_strain5.add_epitope(AG2_ep2_strain5); V5_strain5.add_epitope(AG2_ep3_strain5); V5_strain5.add_epitope(AG2_ep4_strain5); V5_strain5.add_epitope(AG2_ep5_strain5); V5_strain5.add_epitope(AG2_ep6_strain5)
#virus 6
V6_strain6 = func.AntigenType("Type17", 17)
V6_strain6.add_epitope(AG1_ep1_strain6); V6_strain6.add_epitope(AG1_ep2_strain6); V6_strain6.add_epitope(AG1_ep3_strain6); V6_strain6.add_epitope(AG1_ep4_strain6); V6_strain6.add_epitope(AG2_ep1_strain6); V6_strain6.add_epitope(AG2_ep2_strain6); V6_strain6.add_epitope(AG2_ep3_strain6); V6_strain6.add_epitope(AG2_ep4_strain6); V6_strain6.add_epitope(AG2_ep5_strain6); V6_strain6.add_epitope(AG2_ep6_strain6)
#virus 7
V7_strain7 = func.AntigenType("Type18", 18)
V7_strain7.add_epitope(AG1_ep1_strain7); V7_strain7.add_epitope(AG1_ep2_strain7); V7_strain7.add_epitope(AG1_ep3_strain7); V7_strain7.add_epitope(AG1_ep4_strain7); V7_strain7.add_epitope(AG2_ep1_strain7); V7_strain7.add_epitope(AG2_ep2_strain7); V7_strain7.add_epitope(AG2_ep3_strain7); V7_strain7.add_epitope(AG2_ep4_strain7); V7_strain7.add_epitope(AG2_ep5_strain7); V7_strain7.add_epitope(AG2_ep6_strain7)
#virus 8
V8_strain8 = func.AntigenType("Type19", 19)
V8_strain8.add_epitope(AG1_ep1_strain8); V8_strain8.add_epitope(AG1_ep2_strain8); V8_strain8.add_epitope(AG1_ep3_strain8); V8_strain8.add_epitope(AG1_ep4_strain8); V8_strain8.add_epitope(AG2_ep1_strain8); V8_strain8.add_epitope(AG2_ep2_strain8); V8_strain8.add_epitope(AG2_ep3_strain8); V8_strain8.add_epitope(AG2_ep4_strain8); V8_strain8.add_epitope(AG2_ep5_strain8); V8_strain8.add_epitope(AG2_ep6_strain8)
#virus 9
V9_strain9 = func.AntigenType("Type20", 20)
V9_strain9.add_epitope(AG1_ep1_strain9); V9_strain9.add_epitope(AG1_ep2_strain9); V9_strain9.add_epitope(AG1_ep3_strain9); V9_strain9.add_epitope(AG1_ep4_strain9); V9_strain9.add_epitope(AG2_ep1_strain9); V9_strain9.add_epitope(AG2_ep2_strain9); V9_strain9.add_epitope(AG2_ep3_strain9); V9_strain9.add_epitope(AG2_ep4_strain9); V9_strain9.add_epitope(AG2_ep5_strain9); V9_strain9.add_epitope(AG2_ep6_strain9)
#virus 10
V10_strain10 = func.AntigenType("Type21", 21)
V10_strain10.add_epitope(AG1_ep1_strain10); V10_strain10.add_epitope(AG1_ep2_strain10); V10_strain10.add_epitope(AG1_ep3_strain10); V10_strain10.add_epitope(AG1_ep4_strain10); V10_strain10.add_epitope(AG2_ep1_strain10); V10_strain10.add_epitope(AG2_ep2_strain10); V10_strain10.add_epitope(AG2_ep3_strain10); V10_strain10.add_epitope(AG2_ep4_strain10); V10_strain10.add_epitope(AG2_ep5_strain10); V10_strain10.add_epitope(AG2_ep6_strain10)
#virus 11
V11_strain11 = func.AntigenType("Type22", 22)
V11_strain11.add_epitope(AG1_ep1_strain11); V11_strain11.add_epitope(AG1_ep2_strain11); V11_strain11.add_epitope(AG1_ep3_strain11); V11_strain11.add_epitope(AG1_ep4_strain11); V11_strain11.add_epitope(AG2_ep1_strain11); V11_strain11.add_epitope(AG2_ep2_strain11); V11_strain11.add_epitope(AG2_ep3_strain11); V11_strain11.add_epitope(AG2_ep4_strain11); V11_strain11.add_epitope(AG2_ep5_strain11); V11_strain11.add_epitope(AG2_ep6_strain11)
#virus 12
V12_strain12 = func.AntigenType("Type23", 23)
V12_strain12.add_epitope(AG1_ep1_strain12); V12_strain12.add_epitope(AG1_ep2_strain12); V12_strain12.add_epitope(AG1_ep3_strain12); V12_strain12.add_epitope(AG1_ep4_strain12); V12_strain12.add_epitope(AG2_ep1_strain12); V12_strain12.add_epitope(AG2_ep2_strain12); V12_strain12.add_epitope(AG2_ep3_strain12); V12_strain12.add_epitope(AG2_ep4_strain12); V12_strain12.add_epitope(AG2_ep5_strain12); V12_strain12.add_epitope(AG2_ep6_strain12)
##################################
#Add Antigens & Viruses to List. #
##################################
antigen_list = []
antigen_list.append(AG1_strain1); antigen_list.append(AG1_strain2); antigen_list.append(AG1_strain3); antigen_list.append(AG1_strain4); antigen_list.append(AG1_strain5); antigen_list.append(AG1_strain6); antigen_list.append(AG1_strain7); antigen_list.append(AG1_strain8); antigen_list.append(AG1_strain9); antigen_list.append(AG1_strain10); antigen_list.append(AG1_strain11); antigen_list.append(AG1_strain12); antigen_list.append(V1_strain1); antigen_list.append(V2_strain2); antigen_list.append(V3_strain3); antigen_list.append(V4_strain4); antigen_list.append(V5_strain5); antigen_list.append(V6_strain6); antigen_list.append(V7_strain7); antigen_list.append(V8_strain8); antigen_list.append(V9_strain9); antigen_list.append(V10_strain10); antigen_list.append(V11_strain11); antigen_list.append(V12_strain12)
################
#Prime System. #
################
#set up antigen population
V0 = func.Antigen("V0", 0, AG1_strain1)
V1 = func.Antigen("V1", 0, AG1_strain2)
V2 = func.Antigen("V2", 0, AG1_strain3)
V3 = func.Antigen("V3", 0, AG1_strain4)
V4 = func.Antigen("V4", 0, AG1_strain5)
V5 = func.Antigen("V5", 0, AG1_strain6)
V6 = func.Antigen("V6", 0, AG1_strain7)
V7 = func.Antigen("V7", 0, AG1_strain8)
V8 = func.Antigen("V8", 0, AG1_strain9)
V9 = func.Antigen("V9", 0, AG1_strain10)
V10 = func.Antigen("V10", 0, AG1_strain11)
V11 = func.Antigen("V11", 0, AG1_strain12)
#set up virus population
F0 = func.Antigen("F0", 0, V1_strain1)
F1 = func.Antigen("F1", 0, V2_strain2)
F2 = func.Antigen("F2", 0, V3_strain3)
F3 = func.Antigen("F3", 0, V4_strain4)
F4 = func.Antigen("F4", 0, V5_strain5)
F5 = func.Antigen("F5", 0, V6_strain6)
F6 = func.Antigen("F6", 0, V7_strain7)
F7 = func.Antigen("F7", 0, V8_strain8)
F8 = func.Antigen("F8", 0, V9_strain9)
F9 = func.Antigen("F9", 0, V10_strain10)
F10 = func.Antigen("F10", 0, V11_strain11)
F11 = func.Antigen("F11", 0, V12_strain12)
#set up B cell populations
nB_gc = func.BCell("GC_B", 0, antigen_list)
nB_stim = func.BCell("Stimulated_B", 0, antigen_list)
nB_me = func.BCell("Memory_B", 0, antigen_list)
nB_pl = func.BCell("Plasma_B", 0, antigen_list)
nAB = func.BCell("Antibody", 0, antigen_list)
nB = func.BCell("Naive_B", params.Bcell_initial, antigen_list)
nB_llpl = func.BCell("LL_Plasma_B", 0, antigen_list)
GC_population = func.GroupPopulation("GC Bcell Population")
GC_population.add_population(nB_gc)
GC_population.add_population(nB_stim)
#set up population list and add all agents in simulation
PopulationList = []
PopulationList.append(nB_gc)#0
PopulationList.append(nB_stim)#1
PopulationList.append(nB_me)#2
PopulationList.append(nB_pl)#3
PopulationList.append(nAB)#4
PopulationList.append(nB)#5
PopulationList.append(nB_llpl)#6
PopulationList.append(V0)#7
PopulationList.append(V1)#8
PopulationList.append(V2)#9
PopulationList.append(V3)#10
PopulationList.append(V4)#11
PopulationList.append(V5)#12
PopulationList.append(V6)#13
PopulationList.append(V7)#14
PopulationList.append(V8)#15
PopulationList.append(V9)#16
PopulationList.append(V10)#17
PopulationList.append(V11)#18
PopulationList.append(F0) #19
PopulationList.append(F1)#20
PopulationList.append(F2)#21
PopulationList.append(F3)#22
PopulationList.append(F4)#23
PopulationList.append(F5)#24
PopulationList.append(F6)#25
PopulationList.append(F7)#26
PopulationList.append(F8)#27
PopulationList.append(F9)#28
PopulationList.append(F10)#29
PopulationList.append(F11)#30
###############################################
#Define Equations (rates, counts, time, etc). #
###############################################
#Description: Circulating B cell antigen stimulation
eq1a = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), V0, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1b = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), V1, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1c = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), V2, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1d = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), V3, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1e = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), V4, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1f = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), V5, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1g = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), V6, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1h = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), V7, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1i = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), V8, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1j = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), V9, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1k = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), V10, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1l = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), V11, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1m = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), F0, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1n = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), F1, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1o = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), F2, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1p = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), F3, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1q = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), F4, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1r = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), F5, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1s = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), F6, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1t = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), F7, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1u = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), F8, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1v = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), F9, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1w = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), F10, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq1x = func.Stimulation("Free Naive B Cell Stimulation", float(params.stimulation / 10.0), F11, nB, nB_gc, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
#Description: Germinal Center B cell antigen stimulation. Binding affinity with antigen scales exponentially with Hamming Distance
eq2a = func.Stimulation("GC B Cell Stimulation", params.stimulation, V0, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")  #max stimulate rate half life 15min (0.25hr)
eq2b = func.Stimulation("GC B Cell Stimulation", params.stimulation, V1, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
eq2c = func.Stimulation("GC B Cell Stimulation", params.stimulation, V2, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")  #max stimulate rate half life 15min (0.25hr)
eq2d = func.Stimulation("GC B Cell Stimulation", params.stimulation, V3, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
eq2e = func.Stimulation("GC B Cell Stimulation", params.stimulation, V4, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")  #max stimulate rate half life 15min (0.25hr)
eq2f = func.Stimulation("GC B Cell Stimulation", params.stimulation, V5, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
eq2g = func.Stimulation("GC B Cell Stimulation", params.stimulation, V6, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")  #max stimulate rate half life 15min (0.25hr)
eq2h = func.Stimulation("GC B Cell Stimulation", params.stimulation, V7, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
eq2i = func.Stimulation("GC B Cell Stimulation", params.stimulation, V8, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")  #max stimulate rate half life 15min (0.25hr)
eq2j = func.Stimulation("GC B Cell Stimulation", params.stimulation, V9, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
eq2k = func.Stimulation("GC B Cell Stimulation", params.stimulation, V10, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")  #max stimulate rate half life 15min (0.25hr)
eq2l = func.Stimulation("GC B Cell Stimulation", params.stimulation, V11, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
eq2m = func.Stimulation("GC B Cell Stimulation", params.stimulation, F0, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
eq2n = func.Stimulation("GC B Cell Stimulation", params.stimulation, F1, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
eq2o = func.Stimulation("GC B Cell Stimulation", params.stimulation, F2, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
eq2p = func.Stimulation("GC B Cell Stimulation", params.stimulation, F3, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
eq2q = func.Stimulation("GC B Cell Stimulation", params.stimulation, F4, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
eq2r = func.Stimulation("GC B Cell Stimulation", params.stimulation, F5, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
eq2s = func.Stimulation("GC B Cell Stimulation", params.stimulation, F6, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
eq2t = func.Stimulation("GC B Cell Stimulation", params.stimulation, F7, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
eq2u = func.Stimulation("GC B Cell Stimulation", params.stimulation, F8, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
eq2v = func.Stimulation("GC B Cell Stimulation", params.stimulation, F9, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
eq2w = func.Stimulation("GC B Cell Stimulation", params.stimulation, F10, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
eq2x = func.Stimulation("GC B Cell Stimulation", params.stimulation, F11, nB_gc, nB_stim, float(params.tau / 0.25), params.Bcell_aff, "immunogenicity")
#Description: Circulating memory cell stimulation rate
eq11a = func.Stimulation("Memory B cell Simulation", params.stimulation, V0, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11b = func.Stimulation("Memory B cell Simulation", params.stimulation, V1, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11c = func.Stimulation("Memory B cell Simulation", params.stimulation, V2, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11d = func.Stimulation("Memory B cell Simulation", params.stimulation, V3, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11e = func.Stimulation("Memory B cell Simulation", params.stimulation, V4, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11f = func.Stimulation("Memory B cell Simulation", params.stimulation, V5, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11g = func.Stimulation("Memory B cell Simulation", params.stimulation, V6, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11h = func.Stimulation("Memory B cell Simulation", params.stimulation, V7, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11i = func.Stimulation("Memory B cell Simulation", params.stimulation, V8, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11j = func.Stimulation("Memory B cell Simulation", params.stimulation, V9, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11k = func.Stimulation("Memory B cell Simulation", params.stimulation, V10, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11l = func.Stimulation("Memory B cell Simulation", params.stimulation, V11, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11m = func.Stimulation("Memory B cell Simulation", params.stimulation, F0, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11n = func.Stimulation("Memory B cell Simulation", params.stimulation, F1, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11o = func.Stimulation("Memory B cell Simulation", params.stimulation, F2, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11p = func.Stimulation("Memory B cell Simulation", params.stimulation, F3, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11q = func.Stimulation("Memory B cell Simulation", params.stimulation, F4, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11r = func.Stimulation("Memory B cell Simulation", params.stimulation, F5, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11s = func.Stimulation("Memory B cell Simulation", params.stimulation, F6, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11t = func.Stimulation("Memory B cell Simulation", params.stimulation, F7, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11u = func.Stimulation("Memory B cell Simulation", params.stimulation, F8, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11v = func.Stimulation("Memory B cell Simulation", params.stimulation, F9, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11w = func.Stimulation("Memory B cell Simulation", params.stimulation, F10, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")
eq11x = func.Stimulation("Memory B cell Simulation", params.stimulation, F11, nB_me, nB_stim, float(params.tau / 24.0), params.Bcell_aff, "immunogenicity")

#Description: Naive B cell formation in bone marrow
n_antigen = float(len(antigen_list))
multiplier = n_antigen - (n_antigen-1.0)*func.AvgCrossreactivity(antigen_list)
formation_rate = 0.216/multiplier
eq3 = func.Formation("Naive B Cell formation", float(params.tau / formation_rate), nB)  #produces 500 cells every 108hrs

#Description: Germinal Center B cell decay, governed by apoptosis. Modeled as a logistic function with a pre-defined GC population carrying capacity
eq4a = func.PopulationDecay("GC B Cell decay", float(params.tau / (params.tau + 1.0)), float(params.tau / 108.0), params.Bcell_carry, GC_population, nB_gc)  #base half life of 4.5 days (108hrs)

#Description: Basal B cell decay rate
eq4b = func.Decay("Naive B cell decay", float(params.tau / 108.0), nB)

#Description: Germinal Center B cell differentiation rate, modeling affinity-dependent T help
eq5 = func.Differentiation("B cell differentiation", 1.0, nB_stim, nB_gc, nB_me, nB_pl, nB_llpl, 1.0)  #cell cycle is params.tau

#Description: Antibody production by Plasma cells
eq6a = func.Production("Antibody Production", 1.0, nB_pl, nAB)
eq6b = func.Production("Antibody Production", 1.0, nB_llpl, nAB)

#Description: Antibody decay rate
eq7 = func.Decay("Antibody Decay", float(params.tau / 360.0), nAB)  #half life of 10 days (360hrs)

#Description: Intrinsic decay and antibody-dependent antigen clearance
eq8a = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), V0)
eq8b = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), V1)
eq8c = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), V2)
eq8d = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), V3)
eq8e = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), V4)
eq8f = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), V5)
eq8g = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), V6)
eq8h = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), V7)
eq8i = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), V8)
eq8j = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), V9)
eq8k = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), V10)
eq8l = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), V11)
eq8m = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), F0)
eq8n = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), F1)
eq8o = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), F2)
eq8p = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), F3)
eq8q = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), F4)
eq8r = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), F5)
eq8s = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), F6)
eq8t = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), F7)
eq8u = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), F8)
eq8v = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), F9)
eq8w = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), F10)
eq8x = func.Decay("Antigen Decay", float(params.tau / params.ag_decay), F11)

#Description: Antibody-dependent antigen clearance
eq9a = func.Clearance("Antigen Clearance", params.clearance, V0, nAB, 10000.0, params.AB_aff, "clearance")
eq9b = func.Clearance("Antigen Clearance", params.clearance, V1, nAB, 10000.0, params.AB_aff, "clearance")
eq9c = func.Clearance("Antigen Clearance", params.clearance, V2, nAB, 10000.0, params.AB_aff, "clearance")
eq9d = func.Clearance("Antigen Clearance", params.clearance, V3, nAB, 10000.0, params.AB_aff, "clearance")
eq9e = func.Clearance("Antigen Clearance", params.clearance, V4, nAB, 10000.0, params.AB_aff, "clearance")
eq9f = func.Clearance("Antigen Clearance", params.clearance, V5, nAB, 10000.0, params.AB_aff, "clearance")
eq9g = func.Clearance("Antigen Clearance", params.clearance, V6, nAB, 10000.0, params.AB_aff, "clearance")
eq9h = func.Clearance("Antigen Clearance", params.clearance, V7, nAB, 10000.0, params.AB_aff, "clearance")
eq9i = func.Clearance("Antigen Clearance", params.clearance, V8, nAB, 10000.0, params.AB_aff, "clearance")
eq9j = func.Clearance("Antigen Clearance", params.clearance, V9, nAB, 10000.0, params.AB_aff, "clearance")
eq9k = func.Clearance("Antigen Clearance", params.clearance, V10, nAB, 10000.0, params.AB_aff, "clearance")
eq9l = func.Clearance("Antigen Clearance", params.clearance, V11, nAB, 10000.0, params.AB_aff, "clearance")
eq9m = func.Clearance("Antigen Clearance", params.clearance, F0, nAB, 10000.0, params.AB_aff, "clearance")
eq9n = func.Clearance("Antigen Clearance", params.clearance, F1, nAB, 10000.0, params.AB_aff, "clearance")
eq9o = func.Clearance("Antigen Clearance", params.clearance, F2, nAB, 10000.0, params.AB_aff, "clearance")
eq9p = func.Clearance("Antigen Clearance", params.clearance, F3, nAB, 10000.0, params.AB_aff, "clearance")
eq9q = func.Clearance("Antigen Clearance", params.clearance, F4, nAB, 10000.0, params.AB_aff, "clearance")
eq9r = func.Clearance("Antigen Clearance", params.clearance, F5, nAB, 10000.0, params.AB_aff, "clearance")
eq9s = func.Clearance("Antigen Clearance", params.clearance, F6, nAB, 10000.0, params.AB_aff, "clearance")
eq9t = func.Clearance("Antigen Clearance", params.clearance, F7, nAB, 10000.0, params.AB_aff, "clearance")
eq9u = func.Clearance("Antigen Clearance", params.clearance, F8, nAB, 10000.0, params.AB_aff, "clearance")
eq9v = func.Clearance("Antigen Clearance", params.clearance, F9, nAB, 10000.0, params.AB_aff, "clearance")
eq9w = func.Clearance("Antigen Clearance", params.clearance, F10, nAB, 10000.0, params.AB_aff, "clearance")
eq9x = func.Clearance("Antigen Clearance", params.clearance, F11, nAB, 10000.0, params.AB_aff, "clearance")

#Description: Plasma cell decay rate
eq10a = func.Decay("Plasma B cell decay", float(params.tau / 72.0), nB_pl)  #half life of 3 days (72hrs)
eq10b = func.Decay("LL Plasma B cell decay", float(params.tau / 4800.0), nB_llpl)  #half life of 200 days (4800hrs)

#Description: Virus Replication
replication_rate = 6 #replicates every 6 hrs
eq12 = func.Replication("Virus Replication", float(params.tau / replication_rate), PopulationList[params.challengeVirus])

#define a system of reactions
A0 = func.TotalReaction()
A0.add_reaction(eq1a); A0.add_reaction(eq1b); A0.add_reaction(eq1c); A0.add_reaction(eq1d); A0.add_reaction(eq1e); A0.add_reaction(eq1f); A0.add_reaction(eq1g); A0.add_reaction(eq1h); A0.add_reaction(eq1i); A0.add_reaction(eq1j); A0.add_reaction(eq1k); A0.add_reaction(eq1l); A0.add_reaction(eq1m); A0.add_reaction(eq1n); A0.add_reaction(eq1o); A0.add_reaction(eq1p); A0.add_reaction(eq1q); A0.add_reaction(eq1r); A0.add_reaction(eq1s); A0.add_reaction(eq1t); A0.add_reaction(eq1u); A0.add_reaction(eq1v); A0.add_reaction(eq1w); A0.add_reaction(eq1x);
A0.add_reaction(eq2a); A0.add_reaction(eq2b); A0.add_reaction(eq2c); A0.add_reaction(eq2d); A0.add_reaction(eq2e); A0.add_reaction(eq2f); A0.add_reaction(eq2g); A0.add_reaction(eq2h); A0.add_reaction(eq2i); A0.add_reaction(eq2j); A0.add_reaction(eq2k); A0.add_reaction(eq2l); A0.add_reaction(eq2m); A0.add_reaction(eq2n); A0.add_reaction(eq2o); A0.add_reaction(eq2p); A0.add_reaction(eq2q); A0.add_reaction(eq2r); A0.add_reaction(eq2s); A0.add_reaction(eq2t); A0.add_reaction(eq2u); A0.add_reaction(eq2v); A0.add_reaction(eq2w); A0.add_reaction(eq2x);
A0.add_reaction(eq3);
A0.add_reaction(eq4a); A0.add_reaction(eq4b);
A0.add_reaction(eq5);
A0.add_reaction(eq6a); A0.add_reaction(eq6b);
A0.add_reaction(eq7);
A0.add_reaction(eq8a); A0.add_reaction(eq8b); A0.add_reaction(eq8c); A0.add_reaction(eq8d); A0.add_reaction(eq8e); A0.add_reaction(eq8f); A0.add_reaction(eq8g); A0.add_reaction(eq8h); A0.add_reaction(eq8i); A0.add_reaction(eq8j); A0.add_reaction(eq8k); A0.add_reaction(eq8l); A0.add_reaction(eq8m); A0.add_reaction(eq8n); A0.add_reaction(eq8o); A0.add_reaction(eq8p); A0.add_reaction(eq8q); A0.add_reaction(eq8r); A0.add_reaction(eq8s); A0.add_reaction(eq8t); A0.add_reaction(eq8u); A0.add_reaction(eq8v); A0.add_reaction(eq8w); A0.add_reaction(eq8x);
A0.add_reaction(eq9a); A0.add_reaction(eq9b); A0.add_reaction(eq9c); A0.add_reaction(eq9d); A0.add_reaction(eq9e); A0.add_reaction(eq9f); A0.add_reaction(eq9g); A0.add_reaction(eq9h); A0.add_reaction(eq9i); A0.add_reaction(eq9j); A0.add_reaction(eq9k); A0.add_reaction(eq9l); A0.add_reaction(eq9m); A0.add_reaction(eq9n); A0.add_reaction(eq9o); A0.add_reaction(eq9p); A0.add_reaction(eq9q); A0.add_reaction(eq9r); A0.add_reaction(eq9s); A0.add_reaction(eq9t); A0.add_reaction(eq9u); A0.add_reaction(eq9v); A0.add_reaction(eq9w); A0.add_reaction(eq9x);
A0.add_reaction(eq10a); A0.add_reaction(eq10b);
A0.add_reaction(eq11a); A0.add_reaction(eq11b); A0.add_reaction(eq11c); A0.add_reaction(eq11d); A0.add_reaction(eq11e); A0.add_reaction(eq11f); A0.add_reaction(eq11g); A0.add_reaction(eq11h); A0.add_reaction(eq11i); A0.add_reaction(eq11j); A0.add_reaction(eq11k); A0.add_reaction(eq11l); A0.add_reaction(eq11m); A0.add_reaction(eq11n); A0.add_reaction(eq11o); A0.add_reaction(eq11p); A0.add_reaction(eq11q); A0.add_reaction(eq11r); A0.add_reaction(eq11s); A0.add_reaction(eq11t); A0.add_reaction(eq11u); A0.add_reaction(eq11v); A0.add_reaction(eq11w); A0.add_reaction(eq11x);
A0.add_reaction(eq12);
######################
#Output File Naming. #
######################
#output file
data_file = 'RSV_Vopt_' + str(random.getrandbits(11)) + '.txt'
print("Writing to data file " + data_file)
#output conditions file
cond_file = data_file + '_conditions' + '.txt'
file1 = open(cond_file,"a")
file1.write(str( "First Dose Interval_" + str(params.vacc_interval1) + "\n"))
file1.write(str( "Second Dose Interval_" + str(params.vacc_interval2) + "\n"))
file1.write(str( "Third Dose Interval_" + str(params.vacc_interval3) + "\n"))
file1.write(str( "Stop Interval_" + str(params.stop_interval) + "\n"))
file1.write(str( "Infection Interval_" + str(params.infect_interval) + "\n"))
file1.write(str( "First Dose Ag Amount_" + str(params.ag1_initialA) + str("_") + str(params.ag1_initialB) + "\n"))
file1.write(str( "Second Dose Ag Amount_" + str(params.ag2_initialA) + str("_") + str(params.ag2_initialB) + "\n"))
file1.write(str( "Third Dose Ag Amount_" + str(params.ag3_initialA) + str("_") + str(params.ag3_initialB) + "\n"))
file1.write(str( "Virus Amount_" + str(params.virus_initial) + "\n"))
file1.write(str( "First Dose Ags_" + str(params.prime_agA) + str("_") + str(params.prime_agB) + "\n"))
file1.write(str( "Second Dose Ags_" + str(params.boost1_agA) + str("_") + str(params.boost1_agB) + "\n"))
file1.write(str( "Third Dose Ags_" + str(params.boost2_agA) + str("_") + str(params.boost2_agB) + "\n"))
file1.write(str( "Virus Ags_" + str(params.challengeVirus) + "\n"))
file1.close()
##################
#Run Simulation. #
##################
output = func.FileOutput(data_file, 0.1, PopulationList, antigen_list)
output.start()

total_time = float(params.vacc_interval1 * 3.0)
t = 0
while ( t <= total_time ):
    dt = A0.MC_TimeStep()
    t += dt
    A0.MC_React()
    output.write(t)

#dose 1
PopulationList[params.prime_agA].increase(params.ag1_initialA); PopulationList[params.prime_agB].increase(params.ag1_initialB)

total_time = total_time + float(params.vacc_interval2 * 3.0)
while ( t <= total_time ):
    dt = A0.MC_TimeStep()
    t += dt
    A0.MC_React()
    output.write(t)

#dose 2
PopulationList[params.boost1_agA].increase(params.ag2_initialA); PopulationList[params.boost1_agB].increase(params.ag2_initialB)

total_time = total_time + float(params.vacc_interval3 * 3.0)
while ( t <= total_time ):
    dt = A0.MC_TimeStep()
    t += dt
    A0.MC_React()
    output.write(t)

#dose 3
PopulationList[params.boost2_agA].increase(params.ag3_initialA); PopulationList[params.boost2_agB].increase(params.ag3_initialB)

total_time = total_time + float(params.infect_interval * 3.0)
while ( t <= total_time ):
	dt = A0.MC_TimeStep()
	t += dt
	A0.MC_React()
	output.write(t)

#infection
if (params.infect == "yes"):
    virus_initial = params.virus_initial
    PopulationList[params.challengeVirus].increase(virus_initial)

#stop
total_time = total_time + float(params.stop_interval * 3.0)
while ( t <= total_time ):
    dt = A0.MC_TimeStep()
    t += dt
    A0.MC_React()
    output.write(t)

output.finish()
