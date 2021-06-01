%% Full Condense Code for Predicting VEI from Global Database 
clear,clc

%Insert Directory for your computer
%Authors Computer directory: '/home/brodsky'
master_direct=['/Users/paulcolosi/Desktop'];

%%%%%% Download Data %%%%%%%%

%Download the eruption data from the data file
eruption_data =readtable('GVP_Eruption_Results_Clean_3.csv');

%Download the volcano data catelog from the data file
volcano_data =readtable('New_GVP_Volcano_List_Holocene.csv');

%% Creating a Baseline for Accuracy and RMSE 

%%%%%%% Initiation Section %%%%%%%%

%Initialize the baseline and STD matrix for varying predictor and VEI threshold:
accuracy_baseline = zeros(4,1,2,2);  
accuracy_STD = zeros(4,1,2,2);

%Initialize Unhistoried baseline and STD matrix varying predictor and VEI threshold:
accuracy_baseline_unhist = zeros(4,1,2,2);  
accuracy_STD_unhist = zeros(4,1,2,2);

%Initialize RMSE baseline matrix varying predictor and VEI threshold:
RMSE_baseline = zeros(4,1,2,2);  

%Initialize unhistoried RMSE baseline matrix varying predictor and VEI threshold:
RMSE_baseline_unhist = zeros(4,1,2,2); 

%Varying VEI thresholds  
for iVEI=0:3
    %Varying Time thresholds
    for itime = 0:1
        
        if itime == 0 
            eruption_data_time = 1500;
        else
            eruption_data_time = min(eruption_data.StartYear);
        end 
        %Changing between VEI thresholds and Removing specific VEI 
        for i_bol = 0:1
            
            if i_bol == 0 
                eruption_data_complete = eruption_data(eruption_data.StartYear >= eruption_data_time,:);
                eruption_data_clean = eruption_data_complete(eruption_data_complete.VEI >= iVEI,:);      
            else 
                eruption_data_complete = eruption_data(eruption_data.StartYear >= eruption_data_time,:);
                eruption_data_clean = eruption_data_complete(eruption_data_complete.VEI ~= iVEI,:);
            end
            
            %Remove NAN values from Eruption_data
            eruption_data_clean_ind = find(isnan(eruption_data_clean.VEI) == 0);
            eruption_data_clean_nan = eruption_data_clean(eruption_data_clean_ind,:);

            %find unique Volcano names and numbers  
            [VolcanoList,~,IVolcanoList]=unique(eruption_data_clean_nan.VolcanoName);
            [VolcanoList_cat,~,IVolcanoList_cat]=unique(volcano_data.VolcanoNumber);

            %total number of unqiue volcano names.
            N=length(VolcanoList);
            
            %Looping to File each eruption with a unqiue Volcano
            for i=1:N

                %finds all of the indexs of each volcano from the eruption_data set 
                J=find(IVolcanoList==i);
                
                %Creates a structure of called Volcano and runs through all of the
                %unqiue names of volcanos and creates a row in a column called name in 
                %the volcano structure and it draws these name from the object
                %volcanoList that was created above from the erutpion _data frame
                Volcano(i).name=VolcanoList(i);
                
                %We are creating a new for loop so that we can create assign a volcano
                %number to the volcano name but we only want this to be done once. so
                %we take the vector J that is all of the indexes that correspond with a
                %volcano name and make j equal that for each iteration of the for loop
                for j=1:length(J),

                    if j==1
                        
                        %This creates a new column in the sturcture volcano called
                        %number and draws that number directly from the eruption data
                        %frame if the j index is 1 and then assigns it to the unqiue
                        %volcano name.
                        Volcano(i).number=eruption_data_clean_nan(J(j),:).VolcanoNumber;

                    end
                    
                    %Row goes into the eruption data frame and takes all of the
                    %eruption events index by J for each loop of j which is the length
                    %of the number of indexes that each volcano has for the number of
                    %eruptions that have occured in the holocene.
                    Row=eruption_data_clean_nan(J(j),:);
                    
                    %Creates a column that in the structure volcano called eruption and
                    %column is made up of structures that contain all the eruption 
                    %events for the unique volcano names. This uses the function
                    %created below so that the structure collects the following
                    %attributes from the eruption data set.
                    Volcano(i).eruption(j)=FileEruption(Row);
                end
            end

            
            %%%%%%%%%%Creating a Baseline for Historied Volcanoes%%%%%%%%%% 
            
            for k=1:length(Volcano)
                
                % Finding Historied Volcanoes
                if (length(Volcano(k).eruption)>1)

                    Volcano(k).currentVEI=Volcano(k).eruption(1).vei;
                    
                else
                    %This makes all of the volcanoes with 1 eruption labeled NaN because
                    %we would not be able to have a test group and then a value to predict
                    Volcano(k).currentVEI=NaN;
                    
                end
            end
            
            %Removing NaN values from Current VEI 
            matrix = [Volcano.currentVEI]';
            indx_matrix = find(isnan(matrix(:,1)) == 0);
            final_matrix = matrix(indx_matrix,:);

            matrix_volcano = Volcano(indx_matrix);

            %Current VEI Vector     
            currentVEI = [matrix_volcano.currentVEI];
            Nvolc=length(currentVEI);
            Nbstrap=500;
            
            %Using SyntheticProb_fun to create synthetic VEI from distribution
            
            %Calculating the Accuracy baseline
            for j_i=1:Nbstrap
                Synth=SyntheticProb_fun(currentVEI,Nvolc);
                Result(j_i)=(sum(Synth==currentVEI))/Nvolc;
            end
            
            %Calculating the RMSE baseline 
            for j_i = 1:Nbstrap

                synth_vei = SyntheticProb_fun(currentVEI,Nvolc);
                RMSE_bl(j_i) = (sum(((synth_vei-currentVEI).^2)/Nvolc)).^(1/2);

            end

            %Accuracy Baseline based on the distribution 
            accuracy_baseline(iVEI+1,:,itime+1,i_bol+1) = mean(Result); 

            %Accuracy Standard Deviation of the Result of the baseline
            accuracy_STD(iVEI+1,:,itime+1,i_bol+1) = std(Result);
            
            %RMSE Baseline based on the distribution
            RMSE_baseline(iVEI+1,:,itime+1,i_bol+1) = mean(RMSE_bl);
            
            

            %%%%%%%%%Creating a Baseline for Unhistoried Volcanoes%%%%%%%%%
            
            for k=1:length(Volcano)
                %Find Unhistoried Volcanoes
                if (length(Volcano(k).eruption)==1)

                    Volcano(k).currentVEI=Volcano(k).eruption(1).vei;
                else
                    %Set Historied Volcanoes to NAN
                    Volcano(k).currentVEI=NaN;
                end
            end
            
            %Removing the NaN values from the current VEI Vector
            matrix_new = [[Volcano.currentVEI]'];
            indx_matrix_new = find(isnan(matrix_new(:,1)) == 0);
            final_matrix_new = matrix_new(indx_matrix_new,:);
            matrix_volcano_new = Volcano(indx_matrix_new);
            
            %Current VEI Vector
            currentVEI = [matrix_volcano_new.currentVEI];
            Nvolc=length(currentVEI);
            Nbstrap=500;
            
            %Creating Synthetic data for baseline from distribution
            
            %Calculating the Accuracy baseline
            for j_i=1:Nbstrap
                Synth=SyntheticProb_fun(currentVEI,Nvolc);
                Result_unhist(j_i)=(sum(Synth==currentVEI))/Nvolc;
            end
            
            %Calculating the RMSE baseline 
            for j_i = 1:Nbstrap

                synth_vei = SyntheticProb_fun(currentVEI,Nvolc);
                RMSE_bl(j_i) = (sum(((synth_vei-currentVEI).^2)/Nvolc)).^(1/2);

            end

            %Baseline based on the distribution 
            accuracy_baseline_unhist(iVEI+1,:,itime+1,i_bol+1) = mean(Result_unhist);

            %Standard Deviation of the Result of the baseline 
            accuracy_STD_unhist(iVEI+1,:,itime+1,i_bol+1) = std(Result_unhist);
            
            %RMSE Baseline based on the distribution
            RMSE_baseline_unhist(iVEI+1,:,itime+1,i_bol+1) = mean(RMSE_bl);
           
        end
        
        clearvars -except eruption_data volcano_data accuracy_baseline_unhist accuracy_STD_unhist accuracy_baseline accuracy_STD RMSE_baseline RMSE_baseline_unhist iVEI master_direct
    end 
end

%Save the baseline from code into the accuracy database for both
%Accuracy and RMSE Baseline 
save(strcat(master_direct,'/VEI_Research/Accuracy_Database/Accuracy_RMSE_baseline','.mat'))

%% Simple Prediction Model for Accuracy and RMSE 

%Initialize precent accuracy matrix varying predictor and VEI threshold:
accuracy_simple_prediction = zeros(4,5,2,2);
Num_volc = zeros(4,5,2,2);

%Initialize RMSE matrix varying predictor and VEI threshold:
RMSE_simple_prediction = zeros(4,5,2,2);

%Varying VEI thresholds
for iVEI=0:3
    
    %Varying time thresholds
    for itime = 0:1
        
        if itime == 0 
            eruption_data_time = 1500;
        else
            eruption_data_time = min(eruption_data.StartYear);
        end
        
        %Changing between VEI thresholds and Removing specific VEI
        for i_bol = 0:1
            if i_bol == 0 
                eruption_data_complete = eruption_data(eruption_data.StartYear >= eruption_data_time,:);
                eruption_data_clean = eruption_data_complete(eruption_data_complete.VEI >= iVEI,:);
        
            else 
                eruption_data_complete = eruption_data(eruption_data.StartYear >= eruption_data_time,:);
                eruption_data_clean = eruption_data_complete(eruption_data_complete.VEI ~= iVEI,:);
            end
            eruption_data_clean_ind = find(isnan(eruption_data_clean.VEI) == 0);
            eruption_data_clean_nan = eruption_data_clean(eruption_data_clean_ind,:);

            %Find Unique Volcano Numbers and Names  
            [VolcanoList,~,IVolcanoList]=unique(eruption_data_clean_nan.VolcanoName);
            [VolcanoList_cat,~,IVolcanoList_cat]=unique(volcano_data.VolcanoNumber);

            %The lenght of the volcano list vector which is the total number of unqiue
            %volcano names.
            N=length(VolcanoList);

            %i is a vector from 1 to the length of the volcano list 
            for i=1:N

                %finds all of the indexs of each volcano from the eruption_data set 
                J=find(IVolcanoList==i);

                Volcano(i).name=VolcanoList(i);

 
                for j=1:length(J),

                    if j==1

                        Volcano(i).number=eruption_data_clean_nan(J(j),:).VolcanoNumber;

                    end
 
                    Row=eruption_data_clean_nan(J(j),:);
 
                    Volcano(i).eruption(j)=FileEruption(Row);
                end
            end

            %%%%%%%%%Simple prediction: Accuracy%%%%%%%%%%
 
            for k=1:length(Volcano),

                %Dind and use Historical Volcanoes 
                if (length(Volcano(k).eruption)>1)

                    %Calculate mean, median, mode, min, max 
                    Predictors(1)= round(mean([Volcano(k).eruption(2:end).vei]));
                    Predictors(2)= round(median([Volcano(k).eruption(2:end).vei]));
                    Predictors(3)= (mode([Volcano(k).eruption(2:end).vei]));
                    Predictors(4)= (min([Volcano(k).eruption(2:end).vei]));
                    Predictors(5)= (max([Volcano(k).eruption(2:end).vei]));
                    
                    Volcano(k).minVEI=min([Volcano(k).eruption(2:end).vei]);
                    Volcano(k).modeVEI=mode([Volcano(k).eruption(2:end).vei]);
                    Volcano(k).maxVEI=max([Volcano(k).eruption(2:end).vei]);
                    Volcano(k).medianVEI=round(median([Volcano(k).eruption(2:end).vei]));
                    Volcano(k).meanVEI=round(mean([Volcano(k).eruption(2:end).vei]));
                    
                    %Save predictors in prediction field
                    Volcano(k).prediction = Predictors;

                    %Most Recent VEI 
                    Volcano(k).currentVEI=Volcano(k).eruption(1).vei;
                else
                    
                    Volcano(k).minVEI=NaN;
                    Volcano(k).modeVEI=NaN;
                    Volcano(k).maxVEI=NaN;
                    Volcano(k).medianVEI=NaN;
                    Volcano(k).meanVEI=NaN;
                     
                    Volcano(k).prediction= NaN(1,5);
                    Volcano(k).currentVEI=NaN;
                end
            end
            
            %Creating the Simple Prediction for Accuracy 
            
            %Initialze results matrix:
            Results = ones(length(Volcano),length(Predictors)+1);

            %This lists the prediction of the current VEI for the current eruption and then
            %list the actually VEI of the current eruption.
            for i_k=1:length(Volcano)

                Results(i_k,1:end-1) = Volcano(i_k).prediction;
                Results(i_k,end) = Volcano(i_k).currentVEI;
            end
            
            %Removing the nan values from the vector 
            indx=find( isnan(Results(:,1)) == 0 & isnan(Results(:,2)) == 0 & isnan(Results(:,3)) == 0 & isnan(Results(:,4)) == 0 & isnan(Results(:,5)) == 0 & isnan(Results(:,6)) == 0);
            Rtmp = Results(indx,:);

            companswer = Rtmp(:,1:end-1)==Rtmp(:,end);

            Totalcorrect = sum(companswer(:,1:end));
            Totaltest = length(companswer);
            
            %Calculate Accuracy of Global Dataset
            Fractioncorrect = Totalcorrect./Totaltest;
            
            %Saving simplete Predicition and Number of Volcanoes correctly
            %predicted 
            accuracy_simple_prediction(iVEI+1,:,itime+1,i_bol+1) = Fractioncorrect;
            Num_volc(iVEI+1,:,itime+1,i_bol+1) = Totaltest;
            
            %%%%%%%%%% Simple prediction: RMSE %%%%%%%%%% 
            
            for iy = 1:5

                if iy == 1
                    VEI = [Volcano.meanVEI]';

                elseif iy == 2
                    VEI = [Volcano.medianVEI]';

                elseif iy == 3
                    VEI = [Volcano.modeVEI]';

                elseif iy == 4
                    VEI = [Volcano.minVEI]';

                else
                    VEI = [Volcano.maxVEI]';

                end
                matrix_predictorVEI = [[Volcano.currentVEI]' VEI];

                indx_matrix = find(isnan(matrix_predictorVEI(:,1)) == 0 & isnan(matrix_predictorVEI(:,2)) == 0);
                final_matrix_predictorVEI = matrix_predictorVEI(indx_matrix,:);

                syn_predictorVEI_matrix = Volcano(indx_matrix);

                while iy == 1
                    syn_predictorVEI = [syn_predictorVEI_matrix.meanVEI]';
                    break
                end

                while iy == 2
                    syn_predictorVEI = [syn_predictorVEI_matrix.medianVEI]';
                    break
                end

                while iy == 3
                    syn_predictorVEI = [syn_predictorVEI_matrix.modeVEI]';
                    break
                end

                while iy == 4
                    syn_predictorVEI = [syn_predictorVEI_matrix.minVEI]';
                    break
                end

                while iy == 5
                    syn_predictorVEI = [syn_predictorVEI_matrix.maxVEI]';
                    break
                end

                syn_predictorVEI_currentVEI = [syn_predictorVEI_matrix.currentVEI]';
                
                %Calulate RMSE
                RMSE_predictorVEI = (sum(((syn_predictorVEI-syn_predictorVEI_currentVEI).^2))/length(syn_predictorVEI_currentVEI)).^(1/2);

                RMSE_simple_prediction(iVEI+1,iy,itime+1,i_bol+1) = RMSE_predictorVEI;
            end

        end
        %Clear all variables expect eruption_data, volcano_data, accuracy_m
        clearvars -except eruption_data volcano_data accuracy_simple_prediction iVEI Num_volc RMSE_simple_prediction master_direct
    end
end
%Saving result for Accarcy from simple predictions 
save(strcat(master_direct,'/VEI_Research/Accuracy_Database/Accuracy_RMSE_simple_prediction','.mat'))

%% Creating Tables of Attributes for Accuracy and RMSE 

% VEI Thresholds 
for iVEI=0:3
    
    for itime = 0:1
        
        %TimeThresholds 
        if itime == 0 
            eruption_data_time = 1500;
        else
            eruption_data_time = min(eruption_data.StartYear);
        end 
        
        %Grouping greater than or equal to seperate from removing VEI 0,1,2,3
        for i_bol = 0:1
            if i_bol == 0 
                eruption_data_complete = eruption_data(eruption_data.StartYear >= eruption_data_time,:);
                eruption_data_clean = eruption_data_complete(eruption_data_complete.VEI >= iVEI,:);
        
            else 
                eruption_data_complete = eruption_data(eruption_data.StartYear >= eruption_data_time,:);
                eruption_data_clean = eruption_data_complete(eruption_data_complete.VEI ~= iVEI,:);
            end     

            eruption_data_clean_ind = find(isnan(eruption_data_clean.VEI) == 0);
            eruption_data_clean_nan = eruption_data_clean(eruption_data_clean_ind,:);
            
            %Use Create_table_attributes_ED Function
            [Tbl_VEI_allcat,Tbl_VEI_cat,Tbl_VEI_dom,Tbl_VEI_etd,Tbl_VEI_ed,Tbl_VEI_morph,Tbl_VEI_num,Tbl_VEI_tec,Tbl_unhistVEI_cat,Tbl_unhistVEI_dom,Tbl_unhistVEI_morph,Tbl_unhistVEI_tec] = Create_table_attributes_ED(eruption_data_clean_nan,volcano_data);

            %Initialize Titles of Tables 
            time_thresh = ["after_1500","all_time"];
            bol_diff = ["threshold","remove"];
            
            %Save Workspace 
            matlab.io.saveVariablesToScript(strcat(master_direct,'/VEI_Research/Table_Attributes_Script/Tablescript_Varibles_',num2str(iVEI),'_',time_thresh(itime+1),'_',bol_diff(i_bol+1),'.m'))

        end
         clearvars -except eruption_data volcano_data iVEI master_direct itime
        
     end
end


%% Train Models for Machine Learning Algorithms

% I will be splitting the Two Section up into Calculations for
% Classification Learners and Regression learners. Note when training retraining the 
% models these two section of the code take a long time to run (several Hours)
%% Compile the Results for the Machine Learning Application For RMSE 
clearvars -except master_direct
clc

%If Retrain is true then code will retrain models
%If Retrain is false then code will output results from trained model

retrainflag = false;

if retrainflag == 1

    %Downloading the Training Model to Access the RMSE 
    %Load a training model for each attribute and then save the type of traing
    %model and the RMSE which will be save in a varible to continues to
    %increase each time 

    table_names_mach_learn = ["allcat","cat","dom","etd","ed","morph","num","tec"];

    Script_names = ["0_after_1500_threshold","0_after_1500_remove","0_all_time_threshold","0_all_time_remove","1_after_1500_threshold","1_after_1500_remove","1_all_time_threshold","1_all_time_remove","2_after_1500_threshold","2_after_1500_remove","2_all_time_threshold","2_all_time_remove","3_after_1500_threshold","3_after_1500_remove","3_all_time_threshold","3_all_time_remove"];

    for j = 1:length(Script_names)

        load(strcat(master_direct,'/VEI_Research/Table_Attributes_Script/Tablescript_Varibles_',Script_names(j),'.mat'))

        %Initialize count variable
        len_i = 200;

        %Hisotried Volcanoe Model and Table of Attributes 
        for k = 1:length(table_names_mach_learn)

            if k == 1
                for i = 1:len_i
                    [model,RMSE] = trainRegressionModel_allcat(Tbl_VEI_allcat);
                    predictability(i) = RMSE;
                    model_var(i) = model;
                end
                %Model Counter 
                [~,I] = min(abs(predictability - mean(predictability)));
                model_var1 = model_var(I);
                model_var_con = model_var1;

            elseif k == 2
                for i = 1:len_i
                    [model,RMSE] = trainRegressionModel_cat(Tbl_VEI_cat);
                    predictability(i) = RMSE;
                    model_var(i) = model;
                end
                %Model Counter 
                [~,I] = min(abs(predictability - mean(predictability)));
                model_var1 = model_var(I);
                model_var_con = cat(1,model_var_con,model_var1);
                
            elseif k == 3
                for i = 1:len_i
                    [model,RMSE] = trainRegressionModel_dom(Tbl_VEI_dom);
                    predictability(i) = RMSE;
                    model_var(i) = model;
                end

                %Model Counter 
                [~,I] = min(abs(predictability - mean(predictability)));
                model_var1 = model_var(I);
                model_var_con = cat(1,model_var_con,model_var1);

            elseif k == 4
                for i = 1:len_i
                    [model,RMSE] = trainRegressionModel_etd(Tbl_VEI_etd);
                    predictability(i) = RMSE;
                    model_var(i) = model;
                end

                %Model Counter 
                [~,I] = min(abs(predictability - mean(predictability)));
                model_var1 = model_var(I);
                model_var_con = cat(1,model_var_con,model_var1);
            
            elseif k == 5
                for i = 1:len_i
                    [model,RMSE] = trainRegressionModel_ed(Tbl_VEI_ed);
                    predictability(i) = RMSE;
                    model_var(i) = model;
                end

                %Model Counter 
                [~,I] = min(abs(predictability - mean(predictability)));
                model_var1 = model_var(I);
                model_var_con = cat(1,model_var_con,model_var1);

                
            elseif k == 6
                for i = 1:len_i
                    [model,RMSE] = trainRegressionModel_morph(Tbl_VEI_morph);
                    predictability(i) = RMSE;
                    model_var(i) = model;
                end

                %Model Counter 
                [~,I] = min(abs(predictability - mean(predictability)));
                model_var1 = model_var(I);
                model_var_con = cat(1,model_var_con,model_var1);
                
            elseif k == 7
                for i = 1:len_i
                    [model,RMSE] = trainRegressionModel_num(Tbl_VEI_num);
                    predictability(i) = RMSE;
                    model_var(i) = model;
                end

                %Model Counter 
                [~,I] = min(abs(predictability - mean(predictability)));
                model_var1 = model_var(I);
                model_var_con = cat(1,model_var_con,model_var1);

            elseif k == 8
                for i = 1:len_i
                    [model,RMSE] = trainRegressionModel_tec(Tbl_VEI_tec);
                    predictability(i) = RMSE;
                    model_var(i) = model;
                end

                %Model Counter 
                [~,I] = min(abs(predictability - mean(predictability)));
                model_var1 = model_var(I);
                model_var_con = cat(1,model_var_con,model_var1);
                
            end
        end

        %Save the Variables from worksapec as a .mat file 
        clearvars -except model_var_con j k Script_names table_names_mach_learn master_direct

        Attributes = {'All Attributes','Categorical Data','Dominant Rock Type','Eruption Time Difference','Eruption Duration','Morphology','All Numerical Data','Tectonic Setting'};

        %Then make these variables above in a sturcture 
        for istruct = 1:length(Attributes)
            RMSE_Model_struct(istruct).Attributes = Attributes(istruct);
            RMSE_Model_struct(istruct).trainedmodel = model_var_con(istruct);
        end

        %Save Structure to file 
        save(strcat(master_direct,'/VEI_Research/RMSE_Database/Models_RMSE/Model_RMSE_',Script_names(j),'.mat'))
    end
    
else

    %%%%% Code below gives results (RMSE) without having to retrain the models %%%%%%%
    
    table_names_mach_learn = ["allcat","cat","dom","etd","ed","morph","num","tec"];

    table_names_mach_learn_unhist = ["cat","dom","morph","tec"];

    Script_names = ["0_after_1500_threshold","0_after_1500_remove","0_all_time_threshold","0_all_time_remove","1_after_1500_threshold","1_after_1500_remove","1_all_time_threshold","1_all_time_remove","2_after_1500_threshold","2_after_1500_remove","2_all_time_threshold","2_all_time_remove","3_after_1500_threshold","3_after_1500_remove","3_all_time_threshold","3_all_time_remove"];

    for j = 1:length(Script_names)

        %Load trained models 
        load(strcat(master_direct,'/VEI_Research/RMSE_Database/Models_RMSE/Model_RMSE_',Script_names(j),'.mat'),'model_var_con')

        %Load Table of Attributes 
        load(strcat(master_direct,'/VEI_Research/Table_Attributes_Script/Tablescript_Varibles_',Script_names(j),'.mat'))

        %Initializing Variable 
        ki = [1:length(table_names_mach_learn)];

        %%%%%%%%% Historied Volcanoes Table of Attributes %%%%%%%%%%% 

        for k = 1:length(ki)
            if k == 1
                yfit = model_var_con(1).predictFcn(Tbl_VEI_allcat(:,2:11));

                indx = find(isnan(yfit) == 0);

                proper_len_currentVEI = yfit(indx);

                currentVEI = [Tbl_VEI_allcat.currentVEI];

                RMSE_predictability = (sum(((yfit-currentVEI).^2),'omitnan')/length(proper_len_currentVEI)).^(1/2);

                RMSE = RMSE_predictability;

                RMSE_var = RMSE;

            elseif k == 2
 
                yfit = model_var_con(2).predictFcn(Tbl_VEI_cat(:,2:4));

                indx = find(isnan(yfit) == 0);

                proper_len_currentVEI = yfit(indx);

                currentVEI = [Tbl_VEI_cat.currentVEI];

                RMSE_predictability = (sum(((yfit-currentVEI).^2),'omitnan')/length(proper_len_currentVEI)).^(1/2);

                RMSE = RMSE_predictability;

                RMSE_var = cat(1,RMSE_var,RMSE);

            elseif k == 3

                yfit = model_var_con(3).predictFcn(Tbl_VEI_dom(:,2));

                indx = find(isnan(yfit) == 0);

                proper_len_currentVEI = yfit(indx);

                currentVEI = [Tbl_VEI_dom.currentVEI];

                RMSE_predictability = (sum(((yfit-currentVEI).^2),'omitnan')/length(proper_len_currentVEI)).^(1/2);

                RMSE = RMSE_predictability;

                RMSE_var = cat(1,RMSE_var,RMSE);

            elseif k == 4

                yfit = model_var_con(4).predictFcn(Tbl_VEI_etd(:,2));

                indx = find(isnan(yfit) == 0);

                proper_len_currentVEI = yfit(indx);

                currentVEI = [Tbl_VEI_etd.currentVEI];

                RMSE_predictability = (sum(((yfit-currentVEI).^2),'omitnan')/length(proper_len_currentVEI)).^(1/2);

                RMSE = RMSE_predictability;

                RMSE_var = cat(1,RMSE_var,RMSE);
                
            elseif k == 5

                yfit = model_var_con(5).predictFcn(Tbl_VEI_ed(:,2));

                indx = find(isnan(yfit) == 0);

                proper_len_currentVEI = yfit(indx);

                currentVEI = [Tbl_VEI_ed.currentVEI];

                RMSE_predictability = (sum(((yfit-currentVEI).^2),'omitnan')/length(proper_len_currentVEI)).^(1/2);

                RMSE = RMSE_predictability;

                RMSE_var = cat(1,RMSE_var,RMSE);

            elseif k == 6

                yfit = model_var_con(6).predictFcn(Tbl_VEI_morph(:,2));

                indx = find(isnan(yfit) == 0);

                proper_len_currentVEI = yfit(indx);

                currentVEI = [Tbl_VEI_morph.currentVEI];

                RMSE_predictability = (sum(((yfit-currentVEI).^2),'omitnan')/length(proper_len_currentVEI)).^(1/2);

                RMSE = RMSE_predictability;

                RMSE_var = cat(1,RMSE_var,RMSE);

            elseif k == 7

                yfit = model_var_con(7).predictFcn(Tbl_VEI_num(:,2:9));

                indx = find(isnan(yfit) == 0);

                proper_len_currentVEI = yfit(indx);

                currentVEI = [Tbl_VEI_num.currentVEI];

                RMSE_predictability = (sum(((yfit-currentVEI).^2),'omitnan')/length(proper_len_currentVEI)).^(1/2);

                RMSE = RMSE_predictability;

                RMSE_var = cat(1,RMSE_var,RMSE);

            else

                yfit = model_var_con(8).predictFcn(Tbl_VEI_tec(:,2));

                indx = find(isnan(yfit) == 0);

                proper_len_currentVEI = yfit(indx);

                currentVEI = [Tbl_VEI_tec.currentVEI];

                RMSE_predictability = (sum(((yfit-currentVEI).^2),'omitnan')/length(proper_len_currentVEI)).^(1/2);

                RMSE = RMSE_predictability;

                RMSE_var = cat(1,RMSE_var,RMSE);

            end
        end
            %%%%%%%%% Unhistoried Volcanoes Table of Attributes %%%%%%%%%%%

        for k = 1:length(table_names_mach_learn_unhist)

            if k == 1

                yfit = model_var_con(2).predictFcn(Tbl_unhistVEI_cat(:,2:4));

                indx = find(isnan(yfit) == 0);

                proper_len_currentVEI = yfit(indx);

                currentVEI = [Tbl_unhistVEI_cat.currentVEI];

                RMSE_predictability = (sum(((yfit-currentVEI).^2),'omitnan')/length(proper_len_currentVEI)).^(1/2);

                RMSE_unhist = RMSE_predictability;

                RMSE_var_unhist = RMSE_unhist;

             elseif k == 2

                yfit = model_var_con(3).predictFcn(Tbl_unhistVEI_dom(:,2));

                indx = find(isnan(yfit) == 0);

                proper_len_currentVEI = yfit(indx);

                currentVEI = [Tbl_unhistVEI_dom.currentVEI];

                RMSE_predictability = (sum(((yfit-currentVEI).^2),'omitnan')/length(proper_len_currentVEI)).^(1/2);

                RMSE_unhist = RMSE_predictability;

                RMSE_var_unhist = cat(1,RMSE_var_unhist,RMSE_unhist);

             elseif k == 3

                yfit = model_var_con(6).predictFcn(Tbl_unhistVEI_morph(:,2));

                indx = find(isnan(yfit) == 0);

                proper_len_currentVEI = yfit(indx);

                currentVEI = [Tbl_unhistVEI_morph.currentVEI];

                RMSE_predictability = (sum(((yfit-currentVEI).^2),'omitnan')/length(proper_len_currentVEI)).^(1/2);

                RMSE_unhist = RMSE_predictability;

                RMSE_var_unhist = cat(1,RMSE_var_unhist,RMSE_unhist);

            else
                
                yfit = model_var_con(8).predictFcn(Tbl_unhistVEI_tec(:,2));

                indx = find(isnan(yfit) == 0);

                proper_len_currentVEI = yfit(indx);

                currentVEI = [Tbl_unhistVEI_tec.currentVEI];

                RMSE_predictability = (sum(((yfit-currentVEI).^2),'omitnan')/length(proper_len_currentVEI)).^(1/2);

                RMSE_unhist = RMSE_predictability;

                RMSE_var_unhist = cat(1,RMSE_var_unhist,RMSE_unhist);

            end
        end
        clearvars -except RMSE_var RMSE_var_unhist master_direct Script_names j k table_names_mach_learn table_names_mach_learn_unhist ki

        %Save  to file 
        save(strcat(master_direct,'/VEI_Research/RMSE_Database/RMSE_',Script_names(j),'.mat'))

    end
end
%% Compile the Results for the Machine Learning Application For Accuracy 
clearvars -except master_direct
clc

%If Retrain is true then code will retrain models
%If Retrain is false then code will output results from trained model

retrainflag = false;

if retrainflag == 1
    
    %Downloading the Training Model to Access the Accuracy 
    %Load a training model for each attribute and then save the type of traing
    %model and the accuracy which will be save in a varible to continues to
    %increase each time 

    table_names_mach_learn = ["allcat","cat","dom","etd","ed","morph","num","tec"];

    Script_names = ["0_after_1500_threshold","0_after_1500_remove","0_all_time_threshold","0_all_time_remove","1_after_1500_threshold","1_after_1500_remove","1_all_time_threshold","1_all_time_remove","2_after_1500_threshold","2_after_1500_remove","2_all_time_threshold","2_all_time_remove","3_after_1500_threshold","3_after_1500_remove","3_all_time_threshold","3_all_time_remove"];

    for j = 1:length(Script_names)

        load(strcat(master_direct,'/VEI_Research/Table_Attributes_Script/Tablescript_Varibles_',Script_names(j),'.mat'))

        %Initialize count variable
        len_i = 200;

        %Hisotried Volcanoe Model and Table of Attributes 
        for k = 1:length(table_names_mach_learn)

            if k == 1
                for i = 1:len_i
                    [model,accuracy] = trainClassifier_allcat(Tbl_VEI_allcat);
                    predictability(i) = accuracy;
                    model_var(i) = model;
                end

                %Model Counter 
                [~,I] = min(abs(predictability - mean(predictability)));
                model_var1 = model_var(I);
                model_var_con = model_var1; 

            elseif k == 2
                for i = 1:len_i
                    [model,accuracy] = trainClassifier_cat(Tbl_VEI_cat);
                    predictability(i) = accuracy;
                    model_var(i) = model;
                end
                
                %Model Counter 
                [~,I] = min(abs(predictability - mean(predictability)));
                model_var1 = model_var(I);
                model_var_con = cat(1,model_var_con,model_var1); 

            elseif k == 3
                for i = 1:len_i
                    [model,accuracy] = trainClassifier_dom(Tbl_VEI_dom);
                    predictability(i) = accuracy;
                    model_var(i) = model;
                end

                %Model Counter
                [~,I] = min(abs(predictability - mean(predictability)));
                model_var1 = model_var(I);
                model_var_con = cat(1,model_var_con,model_var1); 

            elseif k == 4
                for i = 1:len_i
                    [model,accuracy] = trainClassifier_etd(Tbl_VEI_etd);
                    predictability(i) = accuracy;
                    model_var(i) = model;
                end

                %Model Counter 
                [~,I] = min(abs(predictability - mean(predictability)));

                model_var1 = model_var(I);

                model_var_con = cat(1,model_var_con,model_var1); 

            elseif k == 5
                for i = 1:len_i
                    [model,accuracy] = trainClassifier_ed(Tbl_VEI_ed);
                    predictability(i) = accuracy;
                    model_var(i) = model;
                end

                %Model Counter 
                [~,I] = min(abs(predictability - mean(predictability)));

                model_var1 = model_var(I);

                model_var_con = cat(1,model_var_con,model_var1); 

            elseif k == 6
                for i = 1:len_i
                    [model,accuracy] = trainClassifier_morph(Tbl_VEI_morph);
                    predictability(i) = accuracy;
                    model_var(i) = model;
                end

                %Model Counter 
                [~,I] = min(abs(predictability - mean(predictability)));

                model_var1 = model_var(I);

                model_var_con = cat(1,model_var_con,model_var1);


            elseif k == 7
                for i = 1:len_i
                    [model,accuracy] = trainClassifier_num(Tbl_VEI_num);
                    predictability(i) = accuracy;
                    model_var(i) = model;
                end

                %Model Counter 
                [~,I] = min(abs(predictability - mean(predictability)));

                model_var1 = model_var(I);

                model_var_con = cat(1,model_var_con,model_var1);


            elseif k == 8
                for i = 1:len_i
                    [model,accuracy] = trainClassifier_tec(Tbl_VEI_tec);
                    predictability(i) = accuracy;
                    model_var(i) = model;
                end

                %Model Counter 
                [~,I] = min(abs(predictability - mean(predictability)));

                model_var1 = model_var(I);

                model_var_con = cat(1,model_var_con,model_var1);

            end

        end

        %Save the Variables from worksapec as a .mat file 
        clearvars -except model_var_con j Script_names table_names_mach_learn Script_names master_direct

        Attributes = {'All Attributes','Categorical Data','Dominant Rock Type','Eruption Time Difference','Eruption Duration','Morphology','All Numerical Data','Tectonic Setting'};

        %Then make these variables above in a sturcture 
        for istruct = 1:length(Attributes)
            Model_struct(istruct).Attributes = Attributes(istruct);
            Model_struct(istruct).trainedmodel = model_var_con(istruct);
        end

        %Save Structure to file
        save(strcat(master_direct,'/VEI_Research/Accuracy_Database/Models_Accuracy/Model_Accuracy_',Script_names(j),'.mat'))
    end

else
    %%%%%% Nonretraining Code %%%%%%%%%%%

    %Downloading the Training Model to Access the Accuracy 
    %Load a training model for each attribute and then save the type of traing
    %model and the accuracy which will be save in a varible to continues to
    %increase each time 
    
    table_names_mach_learn = ["allcat","cat","dom","etd","ed","morph","num","tec"];
    
    table_names_mach_learn_unhist = ["cat","dom","morph","tec"];
    
    Script_names = ["0_after_1500_threshold","0_after_1500_remove","0_all_time_threshold","0_all_time_remove","1_after_1500_threshold","1_after_1500_remove","1_all_time_threshold","1_all_time_remove","2_after_1500_threshold","2_after_1500_remove","2_all_time_threshold","2_all_time_remove","3_after_1500_threshold","3_after_1500_remove","3_all_time_threshold","3_all_time_remove"];

    for j = 1:length(Script_names)

        %Load trained models 
        load(strcat(master_direct,'/VEI_Research/Accuracy_Database/Models_Accuracy/Model_Accuracy_',Script_names(j),'.mat'),'model_var_con')

        %Load Table of Attributes
        load(strcat(master_direct,'/VEI_Research/Table_Attributes_Script/Tablescript_Varibles_',Script_names(j),'.mat'))

        for k = 1:length(table_names_mach_learn)
            if k == 1

                % Historied Volcanoes
                yfit = model_var_con(1).predictFcn(Tbl_VEI_allcat(:,2:11));

                bol_result = [yfit == Tbl_VEI_allcat.currentVEI];

                sum_bol = sum(bol_result);

                predictability = sum_bol/length(yfit);

                %Predictability
                accuracy = predictability;

                Accuracy_var = accuracy;

                %Number of Volcanoes 
                num_volc = height(Tbl_VEI_allcat);

                Num_Volc_var = num_volc;


            elseif k == 2

                % Historied Volcanoes
                yfit = model_var_con(2).predictFcn(Tbl_VEI_cat(:,2:4));

                bol_result = [yfit == Tbl_VEI_cat.currentVEI];

                sum_bol = sum(bol_result);

                predictability = sum_bol/length(yfit);

                %Predictibility 
                accuracy = predictability;

                Accuracy_var = cat(1,Accuracy_var,accuracy);

                %Number of Volcanoes
                num_volc = height(Tbl_VEI_cat);

                Num_Volc_var = cat(1,Num_Volc_var,num_volc);

            elseif k == 3

                % Historied Volcanoes
                yfit = model_var_con(3).predictFcn(Tbl_VEI_dom(:,2));

                bol_result = [yfit == Tbl_VEI_dom.currentVEI];

                sum_bol = sum(bol_result);

                predictability = sum_bol/length(yfit);

                %Predictibility
                accuracy = predictability;

                Accuracy_var = cat(1,Accuracy_var,accuracy);

                %Number of Volcanoes
                num_volc = height(Tbl_VEI_dom);

                Num_Volc_var = cat(1,Num_Volc_var,num_volc);
                
                %Save the Intrinsic Attributes 
                %Fix this line so that the one function is removed
                All_val = Tbl_VEI_dom(logical(ones(1,num_volc)),{'rocktype'});

                Corr_pred = Tbl_VEI_dom(bol_result,{'rocktype'});
                
                dom_att_all = All_val{:,:};

                dom_att_corr = Corr_pred{:,:};

            elseif k == 4

                % Historied Volcanoes
                yfit = model_var_con(4).predictFcn(Tbl_VEI_etd(:,2));

                bol_result = [yfit == Tbl_VEI_etd.currentVEI];

                sum_bol = sum(bol_result);

                predictability = sum_bol/length(yfit);

                %Predictibility
                accuracy = predictability;

                Accuracy_var = cat(1,Accuracy_var,accuracy);

               %Number of Volcanoes
                num_volc = height(Tbl_VEI_etd);

                Num_Volc_var = cat(1,Num_Volc_var,num_volc);
            
             elseif k == 5

                % Historied Volcanoes
                yfit = model_var_con(5).predictFcn(Tbl_VEI_ed(:,2));

                bol_result = [yfit == Tbl_VEI_ed.currentVEI];

                sum_bol = sum(bol_result);

                predictability = sum_bol/length(yfit);

                %Predictibility
                accuracy = predictability;

                Accuracy_var = cat(1,Accuracy_var,accuracy);

               %Number of Volcanoes
                num_volc = height(Tbl_VEI_etd);

                Num_Volc_var = cat(1,Num_Volc_var,num_volc);

                
            elseif k == 6

                % Historied Volcanoes
                yfit = model_var_con(6).predictFcn(Tbl_VEI_morph(:,2));

                bol_result = [yfit == Tbl_VEI_morph.currentVEI];

                sum_bol = sum(bol_result);

                predictability = sum_bol/length(yfit);

                %Predictibility
                accuracy = predictability;

                Accuracy_var = cat(1,Accuracy_var,accuracy);

                %Number of Volcanoes
                num_volc = height(Tbl_VEI_morph);

                Num_Volc_var = cat(1,Num_Volc_var,num_volc);
                
                %Save the Intrinsic Attributes 
                %Fix this line so that the one function is removed
                All_val = Tbl_VEI_morph(logical(ones(1,num_volc)),{'morphology'});

                Corr_pred = Tbl_VEI_morph(bol_result,{'morphology'});
                
                morph_att_all = All_val{:,:};

                morph_att_corr = Corr_pred{:,:};

            elseif k == 7

                % Historied Volcanoes
                yfit = model_var_con(7).predictFcn(Tbl_VEI_num(:,2:9));

                bol_result = [yfit == Tbl_VEI_num.currentVEI];

                sum_bol = sum(bol_result);

                predictability = sum_bol/length(yfit);

                %Predictibility
                accuracy = predictability;

                Accuracy_var = cat(1,Accuracy_var,accuracy);

                %Number of Volcanoes
                num_volc = height(Tbl_VEI_num);

                Num_Volc_var = cat(1,Num_Volc_var,num_volc);

            elseif k == 8

                % Historied Volcanoes
                yfit = model_var_con(8).predictFcn(Tbl_VEI_tec(:,2));

                bol_result = [yfit == Tbl_VEI_tec.currentVEI];

                sum_bol = sum(bol_result);

                predictability = sum_bol/length(yfit);

                %Predictibility 
                accuracy = predictability;

                Accuracy_var = cat(1,Accuracy_var,accuracy);

                %Number of Volcanoes
                num_volc = height(Tbl_VEI_tec);

                Num_Volc_var = cat(1,Num_Volc_var,num_volc);
                
                %Save the Intrinsic Attributes 
                %Fix this line so that the one function is removed
                All_val = Tbl_VEI_tec(logical(ones(1,num_volc)),{'tectonicsetting'});

                Corr_pred = Tbl_VEI_tec(bol_result,{'tectonicsetting'});
                
                tec_att_all = All_val{:,:};

                tec_att_corr = Corr_pred{:,:};
                
            end  

        end
       %%%%%%%%%%%% For Unhistoried Volcanoes Table of Attributes%%%%%%%%%%

        for k = 1:length(table_names_mach_learn_unhist)

            if k == 1

                yfit = model_var_con(2).predictFcn(Tbl_unhistVEI_cat(:,2:4));
                bol_result = [yfit == Tbl_unhistVEI_cat.currentVEI];

                sum_bol = sum(bol_result);

                predictability = sum_bol/length(yfit);

                %Accuray and Number of Volcaones
                accuracy_unhist = predictability;

                num_volc_unhist = height(Tbl_unhistVEI_cat);

                Accuracy_var_unhist = accuracy_unhist;

                Num_Volc_var_unhist = num_volc_unhist;

             elseif k == 2

                yfit = model_var_con(3).predictFcn(Tbl_unhistVEI_dom(:,2));
                bol_result = [yfit == Tbl_unhistVEI_dom.currentVEI];

                sum_bol = sum(bol_result);

                predictability = sum_bol/length(yfit);

                %Accuray and Number of Volcaones
                accuracy_unhist = predictability;

                num_volc_unhist = height(Tbl_unhistVEI_dom);

                Accuracy_var_unhist = cat(1,Accuracy_var_unhist,accuracy_unhist);

                Num_Volc_var_unhist = cat(1,Num_Volc_var_unhist,num_volc_unhist);

             elseif k == 3

                yfit = model_var_con(6).predictFcn(Tbl_unhistVEI_morph(:,2));
                bol_result = [yfit == Tbl_unhistVEI_morph.currentVEI];

                sum_bol = sum(bol_result);

                predictability = sum_bol/length(yfit);

                %Accuray and Number of Volcaones
                accuracy_unhist = predictability;

                num_volc_unhist = height(Tbl_unhistVEI_morph);

                Accuracy_var_unhist = cat(1,Accuracy_var_unhist,accuracy_unhist);

                Num_Volc_var_unhist = cat(1,Num_Volc_var_unhist,num_volc_unhist);

             else

                yfit = model_var_con(8).predictFcn(Tbl_unhistVEI_tec(:,2));
                bol_result = [yfit == Tbl_unhistVEI_tec.currentVEI];

                sum_bol = sum(bol_result);

                predictability = sum_bol/length(yfit);

                %Accuray and Number of Volcaones
                accuracy_unhist = predictability;

                num_volc_unhist = height(Tbl_unhistVEI_tec);

                Accuracy_var_unhist = cat(1,Accuracy_var_unhist,accuracy_unhist);

                Num_Volc_var_unhist = cat(1,Num_Volc_var_unhist,num_volc_unhist);

            end
        end
        clearvars -except Accuracy_var Accuracy_var_unhist Num_Volc_var Num_Volc_var_unhist intrinsic_att_corr intrinsic_att_all num_att_corr num_att_all morph_att_corr morph_att_all etd_att_corr etd_att_all numerical_att_corr numerical_att_all tec_att_corr tec_att_all dom_att_corr dom_att_all master_direct Script_names j k table_names_mach_learn table_names_mach_learn_unhist

        %Save  to file 
        save(strcat(master_direct,'/VEI_Research/Accuracy_Database/Accuracy_',Script_names(j),'.mat'))        

    end
end
%% Creating the Figures of Accuracy and RMSE Gain

%I have broken these two up into Accuracy and RMSE respectivily with the 
%unhistoried and historied volcanoes following each Category of Accuracy and 
%RMSE 

%% Ploting the difference between baselines and the Accuracy for Historied Volcanoes 
clearvars -except master_direct
close all
clc

%Load the Accuracy from the simple predictions file 
load(strcat(master_direct,'/VEI_Research/Accuracy_Database/Accuracy_RMSE_simple_prediction','.mat'))

%load the baseline from the Accuracy database file 
load(strcat(master_direct,'/VEI_Research/Accuracy_Database/Accuracy_RMSE_baseline','.mat'))

%Concatinating Varibale names with Machine learnign and simple prediction
%names
Accuracy_name_concat = categorical({'All Attributes','Intrinsic Properties of Volcano','Dominant Petrology','Repose Time','Eruption Duration','Morphology','All Numerical Attributes','Tectonic Setting','Mean VEI','Median VEI','Mode VEI','Minimum VEI','Maximum VEI'});

%Reorder the Names for the plot
Accuracy_name_concat_F = reordercats(Accuracy_name_concat,{'All Attributes','Intrinsic Properties of Volcano','Dominant Petrology','Repose Time','Eruption Duration','Morphology','All Numerical Attributes','Tectonic Setting','Mean VEI','Median VEI','Mode VEI','Minimum VEI','Maximum VEI'});

%Initializing Variable
Accuracy = zeros(4,length(Accuracy_name_concat),2);

Accuracy_num = zeros(4,length(Accuracy_name_concat),2);

Accuracy_baseline_diff = zeros(4,length(Accuracy_name_concat),2);

%Initialize counter
nc = 1;

%Loop through VEI thresholds 
for iVEI = 1:4

    %Loop through time thresholds 
    for itime = 1:2
        
        %Values in variables to be plotted 
        accuracy_simple_prediction_F = accuracy_simple_prediction(iVEI,:,itime,1);
        num_simple_prediction_F = Num_volc(iVEI,:,itime,1);
        
        accuracy_baseline_F = accuracy_baseline(iVEI,:,itime,1)*100;

        % Names of the files saved under Accuracy Data
        Script_names = ["0_after_1500_threshold","0_all_time_threshold","1_after_1500_threshold","1_all_time_threshold","2_after_1500_threshold","2_all_time_threshold","3_after_1500_threshold","3_all_time_threshold"];

        %Load the machine learning accuracy from file
        load(strcat(master_direct,'/VEI_Research/Accuracy_Database/Accuracy_',Script_names(nc),'.mat'))
        
        %Machine learning Variable
        accuracy_machine_learning = Accuracy_var';
        
        num_machine_learning = Num_Volc_var';
        
        %Append the variables together 
        Accuracy_concat = 100*[accuracy_machine_learning,accuracy_simple_prediction_F];
        Num_concat = [num_machine_learning,num_simple_prediction_F];

        %Saving Variable
        Accuracy_baseline_diff(iVEI,:,itime) = Accuracy_concat - accuracy_baseline_F;
        
        Accuracy_num(iVEI,:,itime) = Num_concat;
        
        Accuracy(iVEI,:,itime) = 100*[accuracy_machine_learning,accuracy_simple_prediction_F];
        
        nc = nc + 1;
        
    end
end

%Initializing Variable
Accuracy_R = zeros(1,length(Accuracy_name_concat),2);

Accuracy_num_R = zeros(1,length(Accuracy_name_concat),2);

Accuracy_baseline_diff_R = zeros(1,length(Accuracy_name_concat),2);

%Initialize counter
nci = 1;

%Loop through time thresholds 
for itime = 1:2

    %Values in variables to be plotted 
    accuracy_simple_prediction_F = accuracy_simple_prediction(3,:,itime,2);
    
    num_simple_prediction_F = Num_volc(3,:,itime,2);
    
    accuracy_baseline_F = accuracy_baseline(3,:,itime,2)*100;

    % Names of the files saved under Accuracy Data
    Script_names = ["2_after_1500_remove","2_all_time_remove"];

    %Load the machine learning accuracy from file
    load(strcat(master_direct,'/VEI_Research/Accuracy_Database/Accuracy_',Script_names(nci),'.mat'))

    %Machine learning Variable
    accuracy_machine_learning = Accuracy_var';
        
    num_machine_learning = Num_Volc_var';

    %Append the variables together 
    Accuracy_concat = 100*[accuracy_machine_learning,accuracy_simple_prediction_F];
    Num_concat = [num_machine_learning,num_simple_prediction_F];
    
    %Saving Variables for concatination
    Accuracy_baseline_diff_R(1,:,itime) = Accuracy_concat - accuracy_baseline_F;
    
    Accuracy_num_R(1,:,itime) = Num_concat;
    
    Accuracy_R(1,:,itime) = 100*[accuracy_machine_learning,accuracy_simple_prediction_F];
    
    nci = nci + 1;

end 

%Concatinate the Accuracy baseline difference arrays 
Accuracy_baseline_diff_final = cat(1,Accuracy_baseline_diff,Accuracy_baseline_diff_R);

Accuracy_num_final = cat(1,Accuracy_num,Accuracy_num_R);

Accuracy_final = cat(1,Accuracy,Accuracy_R);

%Plotting Figures

figure(1)
%Model Names 
table_name_final = Accuracy_name_concat_F;
%Preset vector of Bar graph
x = [1:length(table_name_final)];

hold on
%Set Font Size 
set(gca,'FontSize',22)
%Plot Bar graph
b = bar(x,Accuracy_baseline_diff_final(:,:,1)');
%Plot vertical line
xl1 = xline(8.5,'--k','Simple Prediction Models','FontSize',18);
xl2 = xline(8.5,'--k','SVM Models','FontSize',18);
hold off

for l = 1:height(Accuracy_num_final) 
    
    xtips1 = b(l).XEndPoints;

    ytips1 = b(l).YEndPoints;

    labels1 = string(Accuracy_num_final(l,:,1));
    
    % Find negative value indices
    Ineg = find(ytips1 < 0);
    Ilogic = ytips1 > 0;

    % If no negative values
    if isempty(Ineg)
        
        text(xtips1,ytips1,labels1,'HorizontalAlignment','right','VerticalAlignment','middle', 'Color', 'w', 'Rotation', 90,'Fontsize',11)
        
    elseif ~isempty(Ineg)
        
        text(xtips1(Ilogic),ytips1(Ilogic),labels1(Ilogic),'HorizontalAlignment','right','VerticalAlignment','middle', 'Color', 'w', 'Rotation', 90,'Fontsize',11)
        text(xtips1(~Ilogic),ytips1(~Ilogic),labels1(~Ilogic),'HorizontalAlignment','left','VerticalAlignment','middle', 'Color', 'w', 'Rotation', 90,'Fontsize',11)
        
    end    
end

%Changing Aspects of Dashed Line 
xl1.LabelOrientation = 'horizontal';
xl2.LabelHorizontalAlignment = 'left';
xl2.LabelOrientation = 'horizontal';

% Change numerical x-labels to categorical data
set(gca,'XTick',[1:1:length(Accuracy_name_concat)])
set(gca,'XTickLabel',string(table_name_final))

%Set Plot features 
title('All VEI Thresholds, After 1500');
ylabel('Accuracy Gain (%)')
ylim([-13,40])
legend('VEI \geq 0','VEI \geq 1','VEI \geq 2','VEI \geq 3','VEI \neq 2','location','northeast')
set(gcf,'units','normalized','outerposition',[0 0 .85 .85])
xtickangle(45)


figure(2)
%Model Names 
table_name_final = Accuracy_name_concat_F;
%Preset vector of Bar graph
x = [1:length(table_name_final)];

hold on
set(gca,'FontSize',22)
%Plot bar graph 
b = bar(x,Accuracy_baseline_diff_final(:,:,2)');
%Plot vertical line
xl1 = xline(8.5,'--k','Simple Prediction Models','FontSize',18);
xl2 = xline(8.5,'--k','SVM Models','FontSize',18);
hold off

for l = 1:height(Accuracy_num_final) 
    
    xtips1 = b(l).XEndPoints;

    ytips1 = b(l).YEndPoints;

    labels1 = string(Accuracy_num_final(l,:,2));
    
    % Find negative value indices
    Ineg = find(ytips1 < 0);
    Ilogic = ytips1 > 0;

    % If no negative values
    if isempty(Ineg)
        
        text(xtips1,ytips1,labels1,'HorizontalAlignment','right','VerticalAlignment','middle', 'Color', 'w', 'Rotation', 90,'Fontsize',10)
        
    elseif ~isempty(Ineg)
        
        text(xtips1(Ilogic),ytips1(Ilogic),labels1(Ilogic),'HorizontalAlignment','right','VerticalAlignment','middle', 'Color', 'w', 'Rotation', 90,'Fontsize',10)
        text(xtips1(~Ilogic),ytips1(~Ilogic),labels1(~Ilogic),'HorizontalAlignment','left','VerticalAlignment','middle', 'Color', 'w', 'Rotation', 90,'Fontsize',10)
        
    end    
end
%Changing Aspects of Dashed Line 
xl1.LabelOrientation = 'horizontal';
xl2.LabelHorizontalAlignment = 'left';
xl2.LabelOrientation = 'horizontal';

% Change numerical x-labels to categorical data
set(gca,'XTick',[1:1:length(Accuracy_name_concat)])
set(gca,'XTickLabel',string(table_name_final))

title('All VEI Thresholds, Holocene')
ylabel('Accuracy Gain (%)')
ylim([-13,42])
legend('VEI \geq 0','VEI \geq 1','VEI \geq 2','VEI \geq 3','VEI \neq 2','location','northeast')
set(gcf,'units','normalized','outerposition',[0 0 .85 .85])
xtickangle(45)

figure(100)
%Model names 
table_name_final = Accuracy_name_concat_F;
%Preset vector of Bar graph
x = [1:length(table_name_final)];

hold on
set(gca,'FontSize',17)
bar(x,Accuracy_final(:,:,1)')
%Plot vertical line
xl1 = xline(8.5,'--k','Simple Prediction Models','FontSize',18);
xl2 = xline(8.5,'--k','SVM Models','FontSize',18);
hold off
%Changing Aspects of Dashed Line 
xl1.LabelOrientation = 'horizontal';
xl2.LabelHorizontalAlignment = 'left';
xl2.LabelOrientation = 'horizontal';

% Change numerical x-labels to categorical data
set(gca,'XTick',[1:1:length(Accuracy_name_concat)])
set(gca,'XTickLabel',string(table_name_final))

title('All VEI Thresholds, After 1500')
ylabel('Total Accuracy (%)')
ylim([0,100])
legend('VEI \geq 0','VEI \geq 1','VEI \geq 2','VEI \geq 3','VEI \neq 2','location','northeast')
set(gcf,'units','normalized','outerposition',[0 0 .85 .85])
xtickangle(45)

figure(101)
%Model names 
table_name_final = Accuracy_name_concat_F;
%Preset vector of Bar graph
x = [1:length(table_name_final)];

hold on
set(gca,'FontSize',17)
bar(x,Accuracy_final(:,:,2)')
%Plot vertical line
xl1 = xline(8.5,'--k','Simple Prediction Models','FontSize',18);
xl2 = xline(8.5,'--k','SVM Models','FontSize',18);
hold off
%Changing Aspects of Dashed Line 
xl1.LabelOrientation = 'horizontal';
xl2.LabelHorizontalAlignment = 'left';
xl2.LabelOrientation = 'horizontal';

% Change numerical x-labels to categorical data
set(gca,'XTick',[1:1:length(Accuracy_name_concat)])
set(gca,'XTickLabel',string(table_name_final))

title('All VEI Thresholds, Holocene')
ylabel('Total Accuracy (%)')
ylim([0,100])
legend('VEI \geq 0','VEI \geq 1','VEI \geq 2','VEI \geq 3','VEI \neq 2','location','northeast')
set(gcf,'units','normalized','outerposition',[0 0 .85 .85])
xtickangle(45)

%% Ploting the difference between baselines and the Accuracy for unhistoried Volcanoes 
clearvars -except master_direct
clc

%load the baseline from the Accuracy database file 
load(strcat(master_direct,'/VEI_Research/Accuracy_Database/Accuracy_RMSE_baseline','.mat'))

%Concatinating Varibale names  
accuracy_machine_learning_name = categorical({'Intrinsic Properties of Volcano','Dominant Petrology','Morphology','Tectonic Setting'});

%Reorder the Names for the plot
Accuracy_name_concat_F = reordercats(accuracy_machine_learning_name,{'Intrinsic Properties of Volcano','Dominant Petrology','Morphology','Tectonic Setting'});

%Initializing Variable
Accuracy = zeros(4,length(accuracy_machine_learning_name),2);

Accuracy_num = zeros(4,length(accuracy_machine_learning_name),2);

Accuracy_baseline_diff = zeros(4,length(accuracy_machine_learning_name),2);

%Initialize counter
nc = 1;

%Loop through VEI thresholds 
for iVEI = 1:4

    %Loop through time thresholds 
    for itime = 1:2

        %Values in variables to be plotted 
        accuracy_baseline_F = accuracy_baseline_unhist(iVEI,:,itime,1)*100;

        % Names of the files saved under Accuracy Data
        Script_names = ["0_after_1500_threshold","0_all_time_threshold","1_after_1500_threshold","1_all_time_threshold","2_after_1500_threshold","2_all_time_threshold","3_after_1500_threshold","3_all_time_threshold"];

        %Load the machine learning accuracy from file
        load(strcat(master_direct,'/VEI_Research/Accuracy_Database/Accuracy_',Script_names(nc),'.mat'))
        
        %Machine learning Variable
        accuracy_machine_learning = Accuracy_var_unhist';
        
        num_machine_learning = Num_Volc_var_unhist';
        
        %Append the variables together 
        Accuracy_concat = 100*[accuracy_machine_learning];
        
        %Saving Variables for concatination
        Accuracy_baseline_diff(iVEI,:,itime) = Accuracy_concat - accuracy_baseline_F;
        
        Accuracy_num(iVEI,:,itime) = num_machine_learning;
        
        Accuracy(iVEI,:,itime) = 100*[accuracy_machine_learning];
        
        nc = nc + 1;
    end
end


%Initializing Variable
Accuracy_R = zeros(1,length(accuracy_machine_learning_name),2);

Accuracy_num_R = zeros(1,length(accuracy_machine_learning_name),2);

Accuracy_baseline_diff_R = zeros(1,length(accuracy_machine_learning_name),2);

%Initialize counter
nci = 1;

%Loop through time thresholds 
for itime = 1:2

    %Values in variables to be plotted 
    accuracy_baseline_F = accuracy_baseline_unhist(3,:,itime,2)*100;

    % Names of the files saved under Accuracy Data
    Script_names = ["2_after_1500_remove","2_all_time_remove"];

    %Load the machine learning accuracy from file
    load(strcat(master_direct,'/VEI_Research/Accuracy_Database/Accuracy_',Script_names(nci),'.mat'))

    %Machine learning Variable
    accuracy_machine_learning = Accuracy_var_unhist';
    num_machine_learning = Num_Volc_var_unhist';
    
    %Append the variables together 
    Accuracy_concat = 100*[accuracy_machine_learning];
    
    %Saving Variables for concatination
    Accuracy_baseline_diff_R(1,:,itime) = Accuracy_concat - accuracy_baseline_F;
    
    Accuracy_num_R(1,:,itime) = num_machine_learning;
    
    Accuracy_R(1,:,itime) = 100*[accuracy_machine_learning];
    
    nci = nci + 1;

end 

%Concatinate the Accuracy baseline difference arrays 
Accuracy_baseline_diff_final = cat(1,Accuracy_baseline_diff,Accuracy_baseline_diff_R);

Accuracy_num_final = cat(1,Accuracy_num,Accuracy_num_R);

Accuracy_final = cat(1,Accuracy,Accuracy_R);

%All categories 
figure(3)
table_name_final = Accuracy_name_concat_F;
hold on
set(gca,'FontSize',20)
b = bar(table_name_final,Accuracy_baseline_diff_final(:,:,1)');
hold off 

for l = 1:height(Accuracy_num_final) 
    
    xtips1 = b(l).XEndPoints;

    ytips1 = b(l).YEndPoints;

    labels1 = string(Accuracy_num_final(l,:,1));
    
    % Find negative value indices
    Ineg = find(ytips1 < 0);
    Ilogic = ytips1 > 0;

    % If no negative values
    if isempty(Ineg)
        
        text(xtips1,ytips1,labels1,'HorizontalAlignment','right','VerticalAlignment','middle', 'Color', 'w', 'Rotation', 90,'Fontsize',14)
        
    elseif ~isempty(Ineg)
        
        text(xtips1(Ilogic),ytips1(Ilogic),labels1(Ilogic),'HorizontalAlignment','right','VerticalAlignment','middle', 'Color', 'w', 'Rotation', 90,'Fontsize',14)
        text(xtips1(~Ilogic),ytips1(~Ilogic),labels1(~Ilogic),'HorizontalAlignment','left','VerticalAlignment','middle', 'Color', 'w', 'Rotation', 90,'Fontsize',14)
        
    end    
end

title('All Unhistoried VEI Thresholds, After 1500')
ylabel('Accuracy Gain (%)')
ylim([0,21])
legend('VEI \geq 0','VEI \geq 1','VEI \geq 2','VEI \geq 3','VEI \neq 2','location','northeast')
set(gcf,'units','normalized','outerposition',[0 0 .85 .85])

figure(4)
table_name_final = Accuracy_name_concat_F;
hold on
set(gca,'FontSize',20)
b = bar(table_name_final,Accuracy_baseline_diff_final(:,:,2)');
hold off 

for l = 1:height(Accuracy_num_final) 
    
    xtips1 = b(l).XEndPoints;

    ytips1 = b(l).YEndPoints;

    labels1 = string(Accuracy_num_final(l,:,2));
    
    % Find negative value indices
    Ineg = find(ytips1 < 0);
    Ilogic = ytips1 > 0;

    % If no negative values
    if isempty(Ineg)
        
        text(xtips1,ytips1,labels1,'HorizontalAlignment','right','VerticalAlignment','middle', 'Color', 'w', 'Rotation', 90,'Fontsize',14)
        
    elseif ~isempty(Ineg)
        
        text(xtips1(Ilogic),ytips1(Ilogic),labels1(Ilogic),'HorizontalAlignment','right','VerticalAlignment','middle', 'Color', 'w', 'Rotation', 90,'Fontsize',14)
        text(xtips1(~Ilogic),ytips1(~Ilogic),labels1(~Ilogic),'HorizontalAlignment','left','VerticalAlignment','middle', 'Color', 'w', 'Rotation', 90,'Fontsize',14)
        
    end    
end

title('All Unhistoried VEI Thresholds, Holocene')
ylabel('Accuracy Gain (%)')
ylim([0,20.5])
legend('VEI \geq 0','VEI \geq 1','VEI \geq 2','VEI \geq 3','VEI \neq 2','location','northeast')
set(gcf,'units','normalized','outerposition',[0 0 .85 .85])

figure(102)
table_name_final = Accuracy_name_concat_F;
hold on
set(gca,'FontSize',20)
bar(table_name_final,Accuracy_final(:,:,1)')
hold off 
title('All Unhistoried VEI Thresholds, After 1500')
ylabel('Total Accuracy (%)')
ylim([0,100])
legend('VEI \geq 0','VEI \geq 1','VEI \geq 2','VEI \geq 3','VEI \neq 2','location','northeast')
set(gcf,'units','normalized','outerposition',[0 0 .85 .85])

figure(103)
table_name_final = Accuracy_name_concat_F;
hold on
set(gca,'FontSize',20)
bar(table_name_final,Accuracy_final(:,:,2)')
hold off 
title('All Unhistoried VEI Thresholds, Holocene')
ylabel('Total Accuracy (%)')
ylim([0,100])
legend('VEI \geq 0','VEI \geq 1','VEI \geq 2','VEI \geq 3','VEI \neq 2','location','northeast')
set(gcf,'units','normalized','outerposition',[0 0 .85 .85])

%% Ploting the difference between baselines and the RMSE for Historied Volcanoes 
clearvars -except master_direct
clc

%Load the RMSE from the simple predictions file 
load(strcat(master_direct,'/VEI_Research/Accuracy_Database/Accuracy_RMSE_simple_prediction','.mat'))

%load the baseline from the RMSE database file 
load(strcat(master_direct,'/VEI_Research/Accuracy_Database/Accuracy_RMSE_baseline','.mat'))

%Concatinating Varibale names with Machine learning and simplete prediction 
RMSE_name_concat = categorical({'All Attributes','Intrinsic Properties of Volcano','Dominant Petrology','Repose Time','Eruption Duration','Morphology','All Numerical Attributes','Tectonic Setting','Mean VEI','Median VEI','Mode VEI','Min VEI','Max VEI'});

%Reorder the Names for the plot
RMSE_name_concat_F = reordercats(RMSE_name_concat,{'All Attributes','Intrinsic Properties of Volcano','Dominant Petrology','Repose Time','Eruption Duration','Morphology','All Numerical Attributes','Tectonic Setting','Mean VEI','Median VEI','Mode VEI','Min VEI','Max VEI'});

%Initializing Variable
RMSE = zeros(4,length(RMSE_name_concat),2);

RMSE_baseline_diff = zeros(4,length(RMSE_name_concat),2);

%Initialize counter
nc = 1;

%Loop through VEI thresholds 
for iVEI = 1:4

    %Loop through time thresholds 
    for itime = 1:2
        
        %Values in variables to be plotted 
        RMSE_simple_prediction_F = RMSE_simple_prediction(iVEI,:,itime,1);
        
        RMSE_baseline_F = RMSE_baseline(iVEI,:,itime,1);

        % Names of the files saved under Accuracy Data
        Script_names = ["0_after_1500_threshold","0_all_time_threshold","1_after_1500_threshold","1_all_time_threshold","2_after_1500_threshold","2_all_time_threshold","3_after_1500_threshold","3_all_time_threshold"];

        %Load the machine learning accuracy from file
        load(strcat(master_direct,'/VEI_Research/RMSE_Database/RMSE_',Script_names(nc),'.mat'))
        
        %Machine learning Variable
        RMSE_machine_learning = RMSE_var';

        %Append the variables together 
        RMSE_concat = [RMSE_machine_learning,RMSE_simple_prediction_F];
        
        RMSE_baseline_diff(iVEI,:,itime) = RMSE_baseline_F - RMSE_concat;
        
        RMSE(iVEI,:,itime) = [RMSE_machine_learning,RMSE_simple_prediction_F];
        
        nc = nc + 1;
        
    end
end

%Initializing Variable
RMSE_R = zeros(1,length(RMSE_name_concat),2);

RMSE_baseline_diff_R = zeros(1,length(RMSE_name_concat),2);

%Initialize counter
nci = 1;

%Loop through time thresholds 
for itime = 1:2

    %Values in variables to be plotted 
    RMSE_simple_prediction_F = RMSE_simple_prediction(3,:,itime,2);

    RMSE_baseline_F = RMSE_baseline(3,:,itime,2);

    % Names of the files saved under Accuracy Data
    Script_names = ["2_after_1500_remove","2_all_time_remove"];

    %Load the machine learning accuracy from file
    load(strcat(master_direct,'/VEI_Research/RMSE_Database/RMSE_',Script_names(nci),'.mat'))

    %Machine learning Variable
    RMSE_machine_learning = RMSE_var';

    %Append the variables together 
    RMSE_concat = [RMSE_machine_learning,RMSE_simple_prediction_F];

    RMSE_baseline_diff_R(1,:,itime) = RMSE_baseline_F - RMSE_concat;
    
    RMSE_R(1,:,itime) = [RMSE_machine_learning,RMSE_simple_prediction_F];
    
    nci = nci + 1;

end 

%Concatinate the RMSE baseline difference arrays 

RMSE_baseline_diff_final = cat(1,RMSE_baseline_diff,RMSE_baseline_diff_R);

RMSE_final = cat(1,RMSE,RMSE_R);

%All Attributes in all VEI thresholds after 1500 
figure(5)
table_name_final = RMSE_name_concat_F;
%Preset vector of Bar graph
x = [1:length(table_name_final)];

hold on
%Set Font Size 
set(gca,'FontSize',22)
%Plot Bar graph
bar(x,RMSE_baseline_diff_final(:,:,1)')
%Plot vertical line
xl1 = xline(8.5,'--k','Simple Prediction Models','FontSize',18);
xl2 = xline(8.5,'--k','SVM Models','FontSize',18);
hold off

%Changing Aspects of Dashed Line 
xl1.LabelOrientation = 'horizontal';
xl2.LabelHorizontalAlignment = 'left';
xl2.LabelOrientation = 'horizontal';

% Change numerical x-labels to categorical data
set(gca,'XTick',[1:1:length(RMSE_name_concat)])
set(gca,'XTickLabel',string(table_name_final))

title('All VEI Thresholds, After 1500')
ylabel('Relative RMSE')
ylim([-0.35,1.15])
legend('VEI \geq 0','VEI \geq 1','VEI \geq 2','VEI \geq 3','VEI \neq 2','location','northeast')
set(gcf,'units','normalized','outerposition',[0 0 .85 .85])
xtickangle(45)

figure(6)
table_name_final = RMSE_name_concat_F;
%Preset vector of Bar graph
x = [1:length(table_name_final)];

hold on
%Set Font Size 
set(gca,'FontSize',22)
%Plot Bar graph
bar(x,RMSE_baseline_diff_final(:,:,2)')
%Plot vertical line
xl1 = xline(8.5,'--k','Simple Prediction Models','FontSize',18);
xl2 = xline(8.5,'--k','SVM Models','FontSize',18);
hold off

%Changing Aspects of Dashed Line 
xl1.LabelOrientation = 'horizontal';
xl2.LabelHorizontalAlignment = 'left';
xl2.LabelOrientation = 'horizontal';

% Change numerical x-labels to categorical data
set(gca,'XTick',[1:1:length(RMSE_name_concat)])
set(gca,'XTickLabel',string(table_name_final))

title('All VEI Thresholds, Holocene')
ylabel('Relative RMSE')
ylim([-0.38,1.4])
legend('VEI \geq 0','VEI \geq 1','VEI \geq 2','VEI \geq 3','VEI \neq 2','location','northeast')
set(gcf,'units','normalized','outerposition',[0 0 .85 .85])
xtickangle(45)

figure(104)
table_name_final = RMSE_name_concat_F;
%Preset vector of Bar graph
x = [1:length(table_name_final)];

hold on
%Set Font Size 
set(gca,'FontSize',22)

%Plot Bar graph
bar(x,RMSE_final(:,:,1)')
%Plot vertical line
xl1 = xline(8.5,'--k','Simple Prediction Models','FontSize',18);
xl2 = xline(8.5,'--k','SVM Models','FontSize',18);
hold off

%Changing Aspects of Dashed Line 
xl1.LabelOrientation = 'horizontal';
xl2.LabelHorizontalAlignment = 'left';
xl2.LabelOrientation = 'horizontal';

% Change numerical x-labels to categorical data
set(gca,'XTick',[1:1:length(RMSE_name_concat)])
set(gca,'XTickLabel',string(table_name_final))
title('All VEI Thresholds, After 1500')
ylabel('RMSE')
ylim([0,2.2])
legend('VEI \geq 0','VEI \geq 1','VEI \geq 2','VEI \geq 3','VEI \neq 2','location','northwest')
set(gcf,'units','normalized','outerposition',[0 0 .85 .85])
xtickangle(45)
 
figure(105)
table_name_final = RMSE_name_concat_F;
%Preset vector of Bar graph
x = [1:length(table_name_final)];

hold on
%Set Font Size 
set(gca,'FontSize',22)

%Plot Bar graph
bar(x,RMSE_final(:,:,2)')
%Plot vertical line
xl1 = xline(8.5,'--k','Simple Prediction Models','FontSize',18);
xl2 = xline(8.5,'--k','SVM Models','FontSize',18);
hold off

%Changing Aspects of Dashed Line 
xl1.LabelOrientation = 'horizontal';
xl2.LabelHorizontalAlignment = 'left';
xl2.LabelOrientation = 'horizontal';

% Change numerical x-labels to categorical data
set(gca,'XTick',[1:1:length(RMSE_name_concat)])
set(gca,'XTickLabel',string(table_name_final))
title('All VEI Thresholds, Holocene')
ylabel('RMSE')
ylim([0,2.8])
legend('VEI \geq 0','VEI \geq 1','VEI \geq 2','VEI \geq 3','VEI \neq 2','location','northwest')
set(gcf,'units','normalized','outerposition',[0 0 .85 .85])
xtickangle(45)

%% Ploting the difference between baselines and the Accuracy for unhistoried Volcanoes 
clearvars -except master_direct
clc

%load the baseline from the Accuracy database file 
load(strcat(master_direct,'/VEI_Research/Accuracy_Database/Accuracy_RMSE_baseline','.mat'))

%Concatinating Varibale names  
RMSE_machine_learning_name = categorical({'Intrinsic Properties of Volcano','Dominant Petrology','Morphology','Tectonic Setting'});

%Reorder the Names for the plot
RMSE_name_concat_F = reordercats(RMSE_machine_learning_name,{'Intrinsic Properties of Volcano','Dominant Petrology','Morphology','Tectonic Setting'});

%Initializing Variable
RMSE = zeros(4,length(RMSE_machine_learning_name),2);

RMSE_baseline_diff = zeros(4,length(RMSE_machine_learning_name),2);

%Initialize counter
nc = 1;

%Loop through VEI thresholds 
for iVEI = 1:4

    %Loop through time thresholds 
    for itime = 1:2

        %Values in variables to be plotted 
        RMSE_baseline_F = RMSE_baseline_unhist(iVEI,:,itime,1);

        % Names of the files saved under Accuracy Data
        Script_names = ["0_after_1500_threshold","0_all_time_threshold","1_after_1500_threshold","1_all_time_threshold","2_after_1500_threshold","2_all_time_threshold","3_after_1500_threshold","3_all_time_threshold"];

        %Load the machine learning accuracy from file
        load(strcat(master_direct,'/VEI_Research/RMSE_Database/RMSE_',Script_names(nc),'.mat'))
        
        %Machine learning Variable
        RMSE_machine_learning = RMSE_var_unhist';
        
        %Saving Variables
        RMSE_baseline_diff(iVEI,:,itime) = RMSE_baseline_F - RMSE_machine_learning;
        
        RMSE(iVEI,:,itime) = RMSE_machine_learning;
        
        nc = nc + 1;
    end
end


%Initializing Variable
RMSE_R = zeros(1,length(RMSE_machine_learning_name),2);

RMSE_baseline_diff_R = zeros(1,length(RMSE_machine_learning_name),2);

%Initialize counter
nci = 1;

%Loop through time thresholds 
for itime = 1:2

    %Values in variables to be plotted 

    RMSE_baseline_F = RMSE_baseline_unhist(3,:,itime,2);

    % Names of the files saved under Accuracy Data
    Script_names = ["2_after_1500_remove","2_all_time_remove"];

    %Load the machine learning accuracy from file
    load(strcat(master_direct,'/VEI_Research/RMSE_Database/RMSE_',Script_names(nci),'.mat'))

    %Machine learning Variable
    RMSE_machine_learning = RMSE_var_unhist';
    
    %Save Variables
    RMSE_baseline_diff_R(1,:,itime) = RMSE_baseline_F - RMSE_machine_learning;
    
    RMSE_R(1,:,itime) = RMSE_machine_learning;

    nci = nci + 1;

end 

%Concatinate the Accuracy baseline difference arrays 

RMSE_baseline_diff_final = cat(1,RMSE_baseline_diff,RMSE_baseline_diff_R);

RMSE_final = cat(1,RMSE,RMSE_R);

%All Tables at All VEI Thresholds after 1500  
figure(7)
table_name_final = RMSE_name_concat_F;
hold on
%Set Font Size 
set(gca,'FontSize',22)
bar(table_name_final,RMSE_baseline_diff_final(:,:,1)')
hold off 
title('All Unhistoried VEI Thresholds, After 1500')
ylabel('Relative RMSE')
ylim([0,0.75])
legend('VEI \geq 0','VEI \geq 1','VEI \geq 2','VEI \geq 3','VEI \neq 2','location','northeast')
set(gcf,'units','normalized','outerposition',[0 0 .85 .85])

%All Tables at All VEI Thresholds all time 
figure(8)
table_name_final = RMSE_name_concat_F;
hold on
%Set Font Size 
set(gca,'FontSize',22)
bar(table_name_final,RMSE_baseline_diff_final(:,:,2)')
hold off
title('All Unhistoried VEI Thresholds, Holocene')
ylabel('Relative RMSE')
ylim([0,1.2])
legend('VEI \geq 0','VEI \geq 1','VEI \geq 2','VEI \geq 3','VEI \neq 2','location','northeast')
set(gcf,'units','normalized','outerposition',[0 0 .85 .85])

figure(106)
table_name_final = RMSE_name_concat_F;
hold on
%Set Font Size 
set(gca,'FontSize',22)
bar(table_name_final,RMSE_final(:,:,1)')
hold off 
title('All Unhistoried VEI Thresholds, After 1500')
ylabel('RMSE')
ylim([0,2.2])
legend('VEI \geq 0','VEI \geq 1','VEI \geq 2','VEI \geq 3','VEI \neq 2','location','northwest')
set(gcf,'units','normalized','outerposition',[0 0 .85 .85])

figure(107)
table_name_final = RMSE_name_concat_F;
hold on
%Set Font Size 
set(gca,'FontSize',22)
bar(table_name_final,RMSE_final(:,:,2)')
hold off 
title('All Unhistoried VEI Thresholds, Holocene')
ylabel('RMSE')
ylim([0,2.5])
legend('VEI \geq 0','VEI \geq 1','VEI \geq 2','VEI \geq 3','VEI \neq 2','location','northwest')
set(gcf,'units','normalized','outerposition',[0 0 .85 .85])

%% Pie Diagram of Morphology 

%Load the machine learning accuracy from file
load(strcat(master_direct,'/VEI_Research/Accuracy_Database/Accuracy_1_all_time_threshold.mat'))

Morph_all = morph_att_all(:,1);

Morph_corr = morph_att_corr(:,1);

%Full distribution
[Morph_List_full,~, Att_indx] = unique(Morph_all);
%Initalize: 
Morph_len_full = zeros(1,6);

for i = 1:length(Morph_List_full)

    Att_find_ind = find(Att_indx == i);

    Morph_len_full(i) = length(Att_find_ind);

end

%Correct distribution
%Initalize: 
Morph_len_corr = zeros(1,6);

[Morph_List_corr,~, Att_indx] = unique(Morph_corr);

for i = 1:length(Morph_List_corr)

    Att_find_ind = find(Att_indx == i);

    Morph_len_corr(i) = length(Att_find_ind);

end 

for istruct = 1:length(Morph_List_full)
    Pie_morph_struct(istruct).Morph_list_full = Morph_List_full(istruct);
    Pie_morph_struct(istruct).Morph_full = Morph_len_full(istruct);
end 

for istruct = 1:length(Morph_List_corr)
    Pie_morph_struct(istruct).Morph_list_corr = Morph_List_corr(istruct);
    Pie_morph_struct(istruct).Morph_corr = Morph_len_corr(istruct);
end


figure(9)
subplot(2,1,1)
X = [Pie_morph_struct.Morph_full];
labels = [Pie_morph_struct.Morph_list_full];
h = pie(X);
title(strcat('Morphology Full Volcano Distribution Holocene, VEI \geq 1'),'FontSize',20)
legend(labels,'location','best','FontSize',20)
set(h(2:2:end),'FontSize',16);

subplot(2,1,2)
X = [Pie_morph_struct.Morph_corr];
labels = [Pie_morph_struct.Morph_list_corr];
h = pie(X);
title(strcat('Morphology Correctly Forecasted Volcano Distribution Holocene, VEI \geq 1'),'FontSize',20)
legend(labels,'location','best','FontSize',20)
set(h(2:2:end),'FontSize',16);
set(gcf,'units','normalized','outerposition',[0 0 1 1])
