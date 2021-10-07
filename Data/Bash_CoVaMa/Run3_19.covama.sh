
#1. run ViReMa.py 

mkdir Run3_19
python2.7 ~/programs/ViReMa_0.7_4CoVaMa/ViReMa.py SIV_R01/SIVB670 ~/programs/SIV-R01/Output/PID_Consensus/Run3_19.consensus.fasta Run3_19/Run3_19 -F --Output_Tag Run3_19 --Output_Dir Run3_19/ --p 8


#2 Run CoVaMa_Make_Matrices.py 
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Make_Matrices.py Run3_19/Run3_19 SIV/AY032751env_padded.fasta  NT --ViReMa_Output1 Run3_19_CoVaMa_Output.txt --Min_Fusion_Coverage 10 --Mode2 Nucs
  
  
#3. Run CoVaMa_Analyse_Matrices.py
python2.7 ~/programs/CoVaMa_0.7/CoVaMa_Analyse_Matrices.py Run3_19/Run3_19.Total_Matrices.py.pi Output/Run3_19.Output.txt Run3_19/ --Min_Coverage 10 --Min_Fusion_Coverage 10 -OutArray -Weighted
