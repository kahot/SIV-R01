
#1. run ViReMa.py 

mkdir Run0_1_
python2.7 ~/programs/ViReMa_0.7_4CoVaMa/ViReMa.py SIV_R01/SIVB670 ~/programs/SIV-R01/Output/PID_Consensus/Run0_1_.consensus.fasta Run0_1_/Run0_1_ -F --Output_Tag Run0_1_ --Output_Dir Run0_1_/ --p 8


#2 Run CoVaMa_Make_Matrices.py 
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Make_Matrices.py Run0_1_/Run0_1_ SIV/AY032751env_padded.fasta  NT --ViReMa_Output1 Run0_1__CoVaMa_Output.txt --Min_Fusion_Coverage 10 --Mode2 Nucs
  
  
#3. Run CoVaMa_Analyse_Matrices.py
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Analyse_Matrices.py Run0_1_/Run0_1_.Total_Matrices.py.pi Output/Run0_1_.Output.txt Run0_1_/ --Min_Coverage 10 --Min_Fusion_Coverage 10 -OutArray -Weighted
