function [ALL_ANAT] = get_Anatomy()
%LOADANATOMY Summary of this function goes here
%   Detailed explanation goes here
BRAIN_CHARS = {'rPCC',[8, -48, 39],'Right Posterior Cingulat Cortext';... % Right Posterior Cingulat Cortext
            'lPCC',[-8, -49, 38],'Left Posterior Cingulat Cortext';... % Left Posterior Cingulat Cortext
            'AMPFC',[0, 47, -2],'Anterior Medial Prefrontal';... % Anterior Medial Prefrontal
            'LIPC',[-45, -67, 37],'Left Inferior Parietal';... % Left Inferior Parietal
            'RIPC',[45, 67, 37],'Right Inferior Preital';... % Right Inferior Preital
            'VMPFC',[6, 30, -9],'Ventro-Medial Prefrontal';... % Ventro-Medial Prefrontal
            'lDPFC',[-39, 34, 37],'Left Dorsolateral Prefrontal (dorsal)';... % Left Dorsolateral Prefrontal (dorsal)
            'rDPFC',[35, 39, 31],'Right Dorsolateral Prefrontal (dorsal)';... % Right Dorsolateral Prefrontal (dorsal)
            'lDPFC',[-46, 38, 8],'Left Dorsolateral Prefrontal (lateral)';... % Left Dorsolateral Prefrontal (lateral)
            'rDPFC',[43, 38, 12],'Right Dorsolateral Prefrontal (lateral)';... % Right Dorsolateral Prefrontal (lateral)
            'rPSM',[41 -27 47],' Right Primary Sensorimotor';... % Right Primary Sensorimotor
            'lPSM',[-40 -27 47],'Left Primary Sensorimotor';... % Left Primary Sensorimotor
            'rPM',[38 -18 45],'Right Primary Motor';... % Right Primary Motor
            'lPM',[-36 -19 48],'Left Primary Motor';... % Left Primary Motor
            'rSA',[15 -33 48],'Rigth Sensory Association';... % Rigth Sensory Association
            'lPM',[-14 -33 48],'Left Sensory Association';... % Left Sensory Association
            'rPMSM',[28 -1 51],'Right Pre-Motor + Supplementary-Motor';... % Right Pre-Motor + Supplementary-Motor
            'lPMSM',[-28 -2 52],'Left Pre-Motor + Supplementary-Motor';... % Left Pre-Motor + Supplementary-Motor
            'rDACC',[6 33 16],'Right Dorsal Anterior Cingulat';... % Right Dorsal Anterior Cingulat
            'lDACC',[-5 39 20],'Left Dorsal Anterior Cingulat';... % Left Dorsal Anterior Cingulat
            }; 
% source. https://bioimagesuiteweb.github.io/webapp/mni2tal.html
% Cite: This application consists of components of the Yale BioImage Suite
% Package. The MNI to Talairach mapping is from Lacadie et al. NeuroImage 
% 2008. The Brodmann area definitions are from the following abstract: 
% C.M. Lacadie, R. K. Fulbright, J. Arora, R.T.Constable, and X. 
% Papademetris. Brodmann Areas defined in MNI space using a new Tracing 
% Tool in BioImage Suite. Human Brain Mapping, 2008. 

            
ALL_ANAT = table(BRAIN_CHARS(:,1),BRAIN_CHARS(:,2),BRAIN_CHARS(:,3),'VariableNames',{'BrainAcronym','coords_MNI','BrainFullName'});
end

