
% Words and pictures assignment 2
% Submitted By:
%	108771077: Sagardeep Mahapatra
%	108800494: Biswaranjan Panda

% inputs: databaseDirectory containing all the images
% outputs: closestmatches -- a model to pass an argument to svmpredict to classify images in testing set
% example usage -- [closestMatches] = get_c_g_red('/Users/tlberg/Desktop/teaching/Fall_12/hw/hw1/images/');

function [closestMatches] = get_c_g_red(databaseDirectory)

training_set = {};
testing_set = {};
training_count = 1;
testing_count = 1;
attributes = {'red','brown','silver','gold','black'};
label = {};
% determine the training set and the testing set here
text_names = [databaseDirectory, 'descr*.txt'];
db_text = dir(text_names);
for i = 1:length(db_text)
    text_path = [databaseDirectory, db_text(i).name];
    file_id = fopen(text_path,'r');
    % collect all the lines in the text file
    db_lines = [];
    while 1
        line = fgetl(file_id);
        if ~ischar(line)
            break;
        end
        db_lines{end+1} = line;
    end 
    fclose(file_id);
    r = db_lines{1};
    % collect all the words of the text     
    db_words = [];
    while ~isempty(r)
        [t,r] = strtok(r);
        db_words{end + 1} = strip_punctuation(t);
    end
    sum = 0;
    lb_indx = 0;
    % compare whether attributes are present in the set of words
    for j = 1:length(attributes)
        for k = 1:length(db_words)
            if strcmpi(db_words{k},attributes{j})
                sum = sum + 1;
                lb_indx = j;
            end            
        end
    end    
    if (sum == 1)
        training_set{training_count} = strrep(strrep(db_text(i).name,'.txt','.jpg'),'descr_','img_');
        label{training_count} = attributes{lb_indx};
        training_count = training_count + 1;
    else
        testing_set{testing_count} = strrep(strrep(db_text(i).name,'.txt','.jpg'),'descr_','img_');
        testing_count = testing_count + 1;        
    end
end

hsv_hist_train = rand(length(training_set),1000);
hsv_hist_test = rand(length(testing_set),1000);
hsv_hist_train = hsv_hist_train .*0;
hsv_hist_test = hsv_hist_test .*0;

% computing hue-saturation-value color histogram for the training set
for k = 1:length(training_set)
    I = imread([databaseDirectory, training_set{k}]);
    I = imresize(I, [32 32]);
    if(size(I,3)==3)
        I = rgb2hsv(im2double(I));
        H = I(:,:,1);
        S = I(:,:,2);
        V = I(:,:,3);

        maxH = max(H(:));
        maxS = max(S(:));
        maxV = max(V(:));

        intvH = maxH/10;
        intvS = maxS/10;
        intvV = maxV/10;

        for i=1:32
            for j=1:32
                I(i,j,1);
                I(i,j,2);
                I(i,j,3);
                h_bin = uint32(I(i,j,1) / intvH) + 1;
                s_bin = uint32(I(i,j,2) / intvS) + 1;
                v_bin = uint32(I(i,j,3) / intvV) + 1;
                if(h_bin > 10)
                    h_bin = 10;
                end
                if(s_bin > 10)
                    s_bin = 10;
                end
                if(v_bin > 10)
                    v_bin = 10;
                end
                hsv_hist_train(k,(((h_bin-1)*100) + ((s_bin-1)*10) + v_bin)) = hsv_hist_train(k,(((h_bin-1)*100) + ((s_bin-1)*10) + v_bin)) + 1;
            end
        end
    end
end

%normalize the training histogram
for k = 1:length(training_set)
    for j = 1:1000
       hsv_hist_train(k,j) = hsv_hist_train(k,j)/1024;
    end
end

% computing hue-saturation-value color histogram for the testing set
for k = 1:length(testing_set)
    I = imread([databaseDirectory, testing_set{k}]);
    I = imresize(I, [32 32]);
    if(size(I,3)==3)
        I = rgb2hsv(im2double(I));
        H = I(:,:,1);
        S = I(:,:,2);
        V = I(:,:,3);

        maxH = max(H(:));
        maxS = max(S(:));
        maxV = max(V(:));

        intvH = maxH/10;
        intvS = maxS/10;
        intvV = maxV/10;

        for i=1:32
            for j=1:32
                I(i,j,1);
                I(i,j,2);
                I(i,j,3);
                h_bin = uint32(I(i,j,1) / intvH) + 1;
                s_bin = uint32(I(i,j,2) / intvS) + 1;
                v_bin = uint32(I(i,j,3) / intvV) + 1;
                if(h_bin > 10)
                    h_bin = 10;
                end
                if(s_bin > 10)
                    s_bin = 10;
                end
                if(v_bin > 10)
                    v_bin = 10;
                end
                hsv_hist_test(k,(((h_bin-1)*100) + ((s_bin-1)*10) + v_bin)) = hsv_hist_test(k,(((h_bin-1)*100) + ((s_bin-1)*10) + v_bin)) + 1;
            end
        end
    end
end

%normalize the testing histogram
for k = 1:length(testing_set)
    for j = 1:1000
       hsv_hist_test(k,j) = hsv_hist_test(k,j)/1024;
    end
end

%divide the training set into 70% training set and 30% tuning set
len_train_set_70 = round(0.7 * length(training_set));
len_tune_set_30 = length(training_set) - len_train_set_70;

%get training label vector for red
blk_trng_lbl_vector = rand(length(training_set));
blk_trng_lbl_vector = blk_trng_lbl_vector .*0;
for k = 1:length(training_set)
    if (strcmp(label{k},'red'))
        blk_trng_lbl_vector(k) = 1;
    else
        blk_trng_lbl_vector(k) = -1;
    end
end

%divide training label vector into training and tuning label vectors
blk_trng_lbl_vect_70 = [];
blk_tune_lbl_vect_30 = [];
hsv_hist_train_70 = rand(len_train_set_70,1000);
hsv_hist_tune_30 = rand(len_tune_set_30,1000);
hsv_hist_train_70 = hsv_hist_train_70 .*0;
hsv_hist_tune_30 = hsv_hist_tune_30 .*0;

for k = 1:len_train_set_70
    blk_trng_lbl_vect_70(k) = blk_trng_lbl_vector(k);
    for l = 1:1000
        hsv_hist_train_70(k,l) = hsv_hist_train(k,l); 
    end
end
for k = 1:len_tune_set_30
    blk_tune_lbl_vect_30(k) = blk_trng_lbl_vector(len_train_set_70+k);
    for l = 1:1000
        hsv_hist_tune_30(k,l) = hsv_hist_train(len_train_set_70+k,l);
    end
end

blk_trng_lbl_vect_70 = blk_trng_lbl_vect_70';
model1 = svmtrain(blk_trng_lbl_vect_70,hsv_hist_train_70,'-c 800 -g 400 -b 1');

blk_tune_lbl_vect_30 = blk_tune_lbl_vect_30';
[predicted_label,accuracy,decision_values] = svmpredict(blk_tune_lbl_vect_30,hsv_hist_tune_30,model1,'-b 1');

closestMatches = model1;




