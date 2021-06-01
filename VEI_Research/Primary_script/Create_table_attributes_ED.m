function [Tbl_VEI_allcat,Tbl_VEI_cat,Tbl_VEI_dom,Tbl_VEI_etd,Tbl_VEI_ed,Tbl_VEI_morph,Tbl_VEI_num,Tbl_VEI_tec,Tbl_unhistVEI_cat,Tbl_unhistVEI_dom,Tbl_unhistVEI_morph,Tbl_unhistVEI_tec] = Create_table_attributes_ED(eruption_data_clean_nan,volcano_data)

%This code produces all of the Tables for each VEI and time threshold

%Inputs: Cleaned Eruption Data and the Volcano Data 
%Outputs: Tables of Attributes for each VEI and Time threshold 
 
    [VolcanoList,~,IVolcanoList]=unique(eruption_data_clean_nan.VolcanoName);
    [VolcanoList_cat,~,~]=unique(volcano_data.VolcanoNumber);

    N=length(VolcanoList);
 
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
        for j=1:length(J)

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

    %This is a for loop that takes the volcano catelog and adds it as a 
    %structure to the volcano data set 
    for l = 1:length(Volcano)

        %finds the index of the volcano with the same volcano number from the
        %Volcanolist_cat vector
        L = find(Volcano(l).number == VolcanoList_cat);
 
        for m = 1:length(L)

            %Rock type, tectonic setting and Morphology 
            Volcano(l).rocktype=volcano_data(L(m),:).DominantRockType;
            Volcano(l).tectonicsetting=volcano_data(L(m),:).TectonicSetting;
            Volcano(l).primaryvolcanotype=volcano_data(L(m),:).PrimaryVolcanoType;
            Volcano(l).region=volcano_data(L(m),:).Region;

        end
    end

    %Find the unique rock types names
    [rocktype,~,rocktype_indx] = unique([Volcano.rocktype]','stable');
    RockType = [rocktype]';

    %Assigning values to the rock types inorder to make them numerical values 
    %This is based on Silica context. 
    rocktype_values = [57.5,46.5,73,63.3,70,NaN,45,56.2,55,44.75,NaN,49];

    %Assign the Silica concentration value to the
    %various rock types 
    for r_num = 1:length(rocktype_values)
       rock_1 = find(rocktype_indx == r_num);
       rocktype_indx(rock_1,:) = rocktype_values(r_num);
    end

    % Make change to the Array of the Volcano Morphologies to consoladate the
    % morphologies

    %Finding the Unqiue Morphology of Volcanoes
    [primaryvolcanotype,~,primaryvolcanotype_indx] = unique([Volcano.primaryvolcanotype]','stable');
    PrimaryVolcanoType = categorical(primaryvolcanotype);

    for q = 1:length(Volcano)
        if Volcano(q).primaryvolcanotype  == "Complex(es)"
            Volcano(q).primaryvolcanotype = 'Stratovolcano';

        elseif Volcano(q).primaryvolcanotype  == "Complex"
            Volcano(q).primaryvolcanotype = 'Stratovolcano';

        elseif Volcano(q).primaryvolcanotype  == "Compound"
            Volcano(q).primaryvolcanotype = 'Stratovolcano';

        elseif Volcano(q).primaryvolcanotype  == "Cone(s)"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif Volcano(q).primaryvolcanotype  == "Stratovolcano(es)"
            Volcano(q).primaryvolcanotype = 'Stratovolcano';

        elseif Volcano(q).primaryvolcanotype  == "Crater rows"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif Volcano(q).primaryvolcanotype  == "Explosion crater(s)"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif  Volcano(q).primaryvolcanotype  == "Fissure vent"
            Volcano(q).primaryvolcanotype = 'Shield';

        elseif Volcano(q).primaryvolcanotype  == "Fissure vent(s)"
            Volcano(q).primaryvolcanotype = 'Shield';

        elseif Volcano(q).primaryvolcanotype  == "Lava cone"
            Volcano(q).primaryvolcanotype = 'Stratovolcano';

        elseif Volcano(q).primaryvolcanotype  == "Lava dome"
            Volcano(q).primaryvolcanotype = 'Stratovolcano';

        elseif Volcano(q).primaryvolcanotype  == "Lava dome(s)"
            Volcano(q).primaryvolcanotype = 'Stratovolcano';

        elseif Volcano(q).primaryvolcanotype  == "Maar"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif Volcano(q).primaryvolcanotype  == "Maar(s)"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif Volcano(q).primaryvolcanotype  == "Pyroclastic cone"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif Volcano(q).primaryvolcanotype  == "Pyroclastic cone(s)"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif Volcano(q).primaryvolcanotype  == "Pyroclastic shield"
            Volcano(q).primaryvolcanotype = 'Caldera';

        elseif Volcano(q).primaryvolcanotype  == "Caldera(s)"
            Volcano(q).primaryvolcanotype = 'Caldera';

        elseif Volcano(q).primaryvolcanotype  == "Shield(s)"
            Volcano(q).primaryvolcanotype = 'Shield';

        elseif Volcano(q).primaryvolcanotype  == "Stratovolcano?"
            Volcano(q).primaryvolcanotype = 'Stratovolcano';

        elseif Volcano(q).primaryvolcanotype  == "Tuff cone(s)"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif Volcano(q).primaryvolcanotype  == "Tuff ring(s)"
            Volcano(q).primaryvolcanotype = 'Field';

        elseif Volcano(q).primaryvolcanotype  == "Volcanic field"
            Volcano(q).primaryvolcanotype = 'Field';

        end
    end

    for k=1:length(Volcano)

        if (length(Volcano(k).eruption)>1)

            %The we are making a column in the volcano structure that tells use
            %the minimum,maximum, median, mean, mode VEI value of a volcano 
            %eruption history from the lastest recorded eruption to the second
            %to most recent eruption

            Volcano(k).minVEI= (min([Volcano(k).eruption(2:end).vei]));
            Volcano(k).modeVEI= (mode([Volcano(k).eruption(2:end).vei]));
            Volcano(k).maxVEI= (max([Volcano(k).eruption(2:end).vei]));
            Volcano(k).medianVEI= round(median([Volcano(k).eruption(2:end).vei]));
            Volcano(k).meanVEI= round(mean([Volcano(k).eruption(2:end).vei]));

            Volcano(k).currentVEI=Volcano(k).eruption(1).vei;

            %Pulling out the values 
            
            Volcano(k).endyear1 = [Volcano(k).eruption(1).endyear];
            Volcano(k).endmonth1 = [Volcano(k).eruption(1).endmonth];
            Volcano(k).endday1 = [Volcano(k).eruption(1).endday];
            
            Volcano(k).startyear2 = [Volcano(k).eruption(2).startyear];
            Volcano(k).startmonth2 = [Volcano(k).eruption(2).startmonth];
            Volcano(k).startday2 = [Volcano(k).eruption(2).startday];
            
            Volcano(k).endyear2 = [Volcano(k).eruption(2).endyear];
            Volcano(k).endmonth2 = [Volcano(k).eruption(2).endmonth];
            Volcano(k).endday2 = [Volcano(k).eruption(2).endday];
            
            %Calculating the difference and Convert all Categories to years 
            %Eruption time elapsed between eruptions 
            Volcano(k).time2 = datenum(Volcano(k).startyear2,Volcano(k).startmonth2,Volcano(k).startday2);
            Volcano(k).time3 = datenum(Volcano(k).endyear1,Volcano(k).endmonth1,Volcano(k).endday1);
            Volcano(k).time4 = datenum(Volcano(k).endyear2,Volcano(k).endmonth2,Volcano(k).endday2);
            
 
            %Eruption Duration 
            Volcano(k).eruptionduration = Volcano(k).time4 - Volcano(k).time2;
            %Eruption time difference 
            Volcano(k).eruptiontimediff = Volcano(k).time3 - Volcano(k).time2;
 
            %The Silica content of the Rock types
            Volcano(k).silica = rocktype_indx(k);
 
 
        else
            %This makes the volcano with a eruption history of length 1 a NaN
            %value for both the minimum and mode columns
 
            %Statisitical Attributes 
            Volcano(k).minVEI=NaN;
            Volcano(k).modeVEI=NaN;
            Volcano(k).maxVEI=NaN;
            Volcano(k).medianVEI=NaN;
            Volcano(k).meanVEI=NaN;
 
            %Calculating Eruption dynamics
            Volcano(k).endyear1 = NaN;
            Volcano(k).startyear2 = NaN;
            Volcano(k).endmonth1 = NaN;
            Volcano(k).startmonth2 = NaN;
            Volcano(k).endday1 = NaN;
            Volcano(k).startday2 = NaN;
 
            Volcano(k).startyear1 = NaN;
            Volcano(k).startmonth1 = NaN;
            Volcano(k).startday1 = NaN;
 
            %Datenum varibles 
            Volcano(k).time1 = NaN;
            Volcano(k).time2 = NaN;
            Volcano(k).time3 = NaN;
 
            %Eruption Dynamic Attributes 
            Volcano(k).eruptiontimediff=NaN;
            Volcano(k).eruptionduration=NaN;
 
            %Rocktype as numerical Variable 
            Volcano(k).silica=NaN;
 
            %The test Variable 
            Volcano(k).currentVEI=NaN;
 
        end
    end
 
    %Removes the [NaN,NaN,NaN] from the datenum varible so the table
    %variable can be concatenated 
    for k = 1:length(Volcano)
        if sum(isnan(Volcano(k).eruptiontimediff))  == 3
 
            Volcano(k).eruptiontimediff = NaN;
 
        end
 
        if sum(isnan(Volcano(k).eruptionduration))  == 3
 
            Volcano(k).eruptionduration = NaN;
 
        end
 
        %Removes one eruption duration error in Data Base 
        if Volcano(k).eruptionduration < 0 
 
            Volcano(k).eruptionduration = NaN;
 
        end    
    end


    time_thresh = ["after_1500","all_time"];

    % Creating Table of Attributes for Machine Learning Application
    for iy = 1:8

        if iy == 1
            matrix_VEI = [[Volcano.currentVEI]' [Volcano.eruptiontimediff]'];
            indx_matrix = find(isnan(matrix_VEI(:,1)) == 0 & isnan(matrix_VEI(:,2)) == 0);
            matrix_volcano = Volcano(indx_matrix);
            %Repose Time 
            Tbl_VEI_etd = table([matrix_volcano.currentVEI]',[matrix_volcano.eruptiontimediff]','VariableNames',{'currentVEI','eruptiontimediff'});
            
        elseif iy == 2
            matrix_VEI = [[Volcano.currentVEI]' [Volcano.eruptionduration]'];
            indx_matrix = find(isnan(matrix_VEI(:,1)) == 0 & isnan(matrix_VEI(:,2)) == 0);
            matrix_volcano = Volcano(indx_matrix);
 
            Tbl_VEI_ed = table([matrix_volcano.currentVEI]',[matrix_volcano.eruptionduration]','VariableNames',{'currentVEI','eruptionduration'});
 
        elseif iy == 3
            matrix_VEI = [[Volcano.currentVEI]' [Volcano.medianVEI]' [Volcano.meanVEI]' [Volcano.modeVEI]' [Volcano.minVEI]' [Volcano.maxVEI]' [Volcano.eruptiontimediff]' [Volcano.eruptionduration]' [Volcano.silica]'];
            indx_matrix = find(isnan(matrix_VEI(:,1)) == 0 & isnan(matrix_VEI(:,2)) == 0 & isnan(matrix_VEI(:,3)) == 0 & isnan(matrix_VEI(:,4)) == 0 & isnan(matrix_VEI(:,5)) == 0 & isnan(matrix_VEI(:,6)) == 0 & isnan(matrix_VEI(:,7)) == 0 & isnan(matrix_VEI(:,8)) == 0 & isnan(matrix_VEI(:,9)) == 0);
            matrix_volcano = Volcano(indx_matrix);
 
            Tbl_VEI_num = table([matrix_volcano.currentVEI]',[matrix_volcano.meanVEI]',[matrix_volcano.medianVEI]',[matrix_volcano.modeVEI]',[matrix_volcano.minVEI]',[matrix_volcano.maxVEI]',[matrix_volcano.eruptiontimediff]',[matrix_volcano.eruptionduration]',[matrix_volcano.silica]','VariableNames',{'currentVEI','meanVEI','medianVEI','modeVEI','minVEI','maxVEI','eruptiontimediff','eruptionduration','silica'});

        elseif iy == 4
            matrix_VEI = [[Volcano.currentVEI]'];
            indx_matrix = find(isnan(matrix_VEI(:,1)) == 0);
            matrix_volcano = Volcano(indx_matrix);
            %Tectonic Setting 
            Tbl_VEI_tec=table([matrix_volcano.currentVEI]',[matrix_volcano.tectonicsetting]','VariableNames',{'currentVEI','tectonicsetting'});

        elseif iy == 5
            matrix_VEI = [[Volcano.currentVEI]'];
            indx1_matrix = find(isnan(matrix_VEI(:,1)) == 0);
            matrix_volcano1 = Volcano(indx1_matrix);

            indx2_matrix = ones(length(matrix_volcano1),1); 
            for q = 1:length(matrix_volcano1)
                if matrix_volcano1(q).rocktype == ""
                    indx2_matrix(q) = 0;
                elseif matrix_volcano1(q).rocktype == "No Data (checked)"
                    indx2_matrix(q) = 0;
                end 
            end
            indx_matrix = find(indx2_matrix == 1);
            matrix_volcano = matrix_volcano1(indx_matrix);
            %Dominant Petrology 
            Tbl_VEI_dom=table([matrix_volcano.currentVEI]',[matrix_volcano.rocktype]','VariableNames',{'currentVEI','rocktype'});


        elseif iy == 6
            matrix_VEI = [[Volcano.currentVEI]'];
            indx_matrix = find(isnan(matrix_VEI(:,1)) == 0);
            matrix_volcano = Volcano(indx_matrix);
            %Morphology 
            Tbl_VEI_morph=table([matrix_volcano.currentVEI]',[matrix_volcano.primaryvolcanotype]','VariableNames',{'currentVEI','morphology'});

        elseif iy == 7
            matrix_VEI = [[Volcano.currentVEI]'];
            indx1_matrix = find(isnan(matrix_VEI(:,1)) == 0);
            matrix_volcano1 = Volcano(indx1_matrix);

            indx2_matrix = ones(length(matrix_volcano1),1); 
            for q = 1:length(matrix_volcano1)
                if matrix_volcano1(q).rocktype == ""
                    indx2_matrix(q) = 0;
                elseif matrix_volcano1(q).rocktype == "No Data (checked)"
                    indx2_matrix(q) = 0;
                end 
            end

            indx_matrix = find(indx2_matrix == 1);
            matrix_volcano = matrix_volcano1(indx_matrix);
            %Intrinisc Properties of Volcanoes 
            Tbl_VEI_cat=table([matrix_volcano.currentVEI]',[matrix_volcano.rocktype]',[matrix_volcano.tectonicsetting]',[matrix_volcano.primaryvolcanotype]','VariableNames',{'currentVEI','rocktype','tectonicsetting','morphology'});    

        else
            matrix_VEI = [[Volcano.currentVEI]' [Volcano.medianVEI]' [Volcano.meanVEI]' [Volcano.modeVEI]' [Volcano.minVEI]' [Volcano.maxVEI]' [Volcano.eruptiontimediff]' [Volcano.eruptionduration]'];
            indx1_matrix = find(isnan(matrix_VEI(:,1)) == 0 & isnan(matrix_VEI(:,2)) == 0 & isnan(matrix_VEI(:,3)) == 0 & isnan(matrix_VEI(:,4)) == 0 & isnan(matrix_VEI(:,5)) == 0 & isnan(matrix_VEI(:,6)) == 0 & isnan(matrix_VEI(:,7)) == 0 & isnan(matrix_VEI(:,8)) == 0);
            matrix_volcano1 = Volcano(indx1_matrix);
 
            indx2_matrix = ones(length(matrix_volcano1),1); 
            for q = 1:length(matrix_volcano1)
                if matrix_volcano1(q).rocktype == ""
                    indx2_matrix(q) = 0;
                elseif matrix_volcano1(q).rocktype == "No Data (checked)"
                    indx2_matrix(q) = 0;
                end 
            end
 
            indx_matrix = find(indx2_matrix == 1);
            matrix_volcano = matrix_volcano1(indx_matrix);
 
            Tbl_VEI_allcat=table([matrix_volcano.currentVEI]',[matrix_volcano.medianVEI]',[matrix_volcano.primaryvolcanotype]',[matrix_volcano.rocktype]',[matrix_volcano.tectonicsetting]',[matrix_volcano.meanVEI]',[matrix_volcano.minVEI]',[matrix_volcano.maxVEI]',[matrix_volcano.modeVEI]',[matrix_volcano.eruptiontimediff]',[matrix_volcano.eruptionduration]','VariableNames',{'currentVEI','medianVEI','morphology','RockType','TectonicSetting','meanVEI','minVEI','maxVEI','modeVEI','eruptiontimediff','eruptionduration'});
        end
    end

    %%%%%%% Making the Unhistoried Table of Attributes %%%%%%% 

    for k=1:length(Volcano)

        if (length(Volcano(k).eruption)==1)
            
            %Current or Most Recent VEI
            Volcano(k).currentVEI=Volcano(k).eruption(1).vei;

        else

            %The test Variable 
            Volcano(k).currentVEI=NaN;

        end
    end

    time_thresh = ["after_1500","all_time"];
    bol_diff = ["threshold","remove"];

    for iy = 1:4

        if iy == 1
            matrix_VEI = [[Volcano.currentVEI]'];
            indx_matrix = find(isnan(matrix_VEI(:,1)) == 0);
            matrix_volcano = Volcano(indx_matrix);
            %Unhistoried Tectonic Setting 
            Tbl_unhistVEI_tec=table([matrix_volcano.currentVEI]',[matrix_volcano.tectonicsetting]','VariableNames',{'currentVEI','tectonicsetting'});

        elseif iy == 2
            matrix_VEI = [[Volcano.currentVEI]'];
            indx1_matrix = find(isnan(matrix_VEI(:,1)) == 0);
            matrix_volcano1 = Volcano(indx1_matrix);

            indx2_matrix = ones(length(matrix_volcano1),1); 
            for q = 1:length(matrix_volcano1)
                if matrix_volcano1(q).rocktype == ""
                    indx2_matrix(q) = 0;
                elseif matrix_volcano1(q).rocktype == "No Data (checked)"
                    indx2_matrix(q) = 0;
                end 
            end

            indx_matrix = find(indx2_matrix == 1);
            matrix_volcano = matrix_volcano1(indx_matrix);
            %Unhistoried Dominant Petrology
            Tbl_unhistVEI_dom=table([matrix_volcano.currentVEI]',[matrix_volcano.rocktype]','VariableNames',{'currentVEI','rocktype'});

        elseif iy == 3
            matrix_VEI = [[Volcano.currentVEI]'];
            indx_matrix = find(isnan(matrix_VEI(:,1)) == 0);
            matrix_volcano = Volcano(indx_matrix);
            %Unhistoried Morphology
            Tbl_unhistVEI_morph=table([matrix_volcano.currentVEI]',[matrix_volcano.primaryvolcanotype]','VariableNames',{'currentVEI','morphology'});

        else
            matrix_VEI = [[Volcano.currentVEI]'];
            indx1_matrix = find(isnan(matrix_VEI(:,1)) == 0);
            matrix_volcano1 = Volcano(indx1_matrix);

            indx2_matrix = ones(length(matrix_volcano1),1); 
            for q = 1:length(matrix_volcano1)
                if matrix_volcano1(q).rocktype == ""
                    indx2_matrix(q) = 0;
                elseif matrix_volcano1(q).rocktype == "No Data (checked)"
                    indx2_matrix(q) = 0;
                end 
            end

            indx_matrix = find(indx2_matrix == 1);
            matrix_volcano = matrix_volcano1(indx_matrix);
            %Unhistoried All Availible Attributes 
            Tbl_unhistVEI_cat=table([matrix_volcano.currentVEI]',[matrix_volcano.rocktype]',[matrix_volcano.tectonicsetting]',[matrix_volcano.primaryvolcanotype]','VariableNames',{'currentVEI','rocktype','tectonicsetting','morphology'});
        end
    end


end

