
% Words and pictures assignment 2
% Submitted By:
%	108771077: Sagardeep Mahapatra
%	108800494: Biswaranjan Panda

% inputs: databaseDirectory containing all the images
% outputs: closestmatches -- a cell array with the filenames of the images forming the training set

% example usage -- [closestMatches] = final_classifier('/Users/tlberg/Desktop/teaching/Fall_12/hw/hw1/images/');

function [closestMatches] = final_classifier(databaseDirectory)

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

%get training label vector for black
for k = 1:length(training_set)
    if (strcmp(label{k},'black'))
        blak_trng_lbl_vector(k) = 1;
    else
        blak_trng_lbl_vector(k) = -1;
    end
end

%get training label vector for red
for k = 1:length(training_set)
    if (strcmp(label{k},'red'))
        red_trng_lbl_vector(k) = 1;
    else
        red_trng_lbl_vector(k) = -1;
    end
end

%get training label vector for brown
for k = 1:length(training_set)
    if (strcmp(label{k},'brown'))
        brown_trng_lbl_vector(k) = 1;
    else
        brown_trng_lbl_vector(k) = -1;
    end
end


%get training label vector for silver
for k = 1:length(training_set)
    if (strcmp(label{k},'silver'))
        silvr_trng_lbl_vector(k) = 1;
    else
        silvr_trng_lbl_vector(k) = -1;
    end
end

%get training label vector for gold
for k = 1:length(training_set)
    if (strcmp(label{k},'gold'))
        gold_trng_lbl_vector(k) = 1;
    else
        gold_trng_lbl_vector(k) = -1;
    end
end


blak_trng_lbl_vector = blak_trng_lbl_vector';
model_blak = svmtrain(blak_trng_lbl_vector,hsv_hist_train,'-c 700 -g 0.01 -b 1');

red_trng_lbl_vector = red_trng_lbl_vector';
model_red = svmtrain(red_trng_lbl_vector,hsv_hist_train,'-c 800 -g 400 -b 1');

silvr_trng_lbl_vector = silvr_trng_lbl_vector';
model_silver = svmtrain(silvr_trng_lbl_vector,hsv_hist_train,'-c 800 -g 600 -b 1');

brown_trng_lbl_vector = brown_trng_lbl_vector';
model_brown = svmtrain(brown_trng_lbl_vector,hsv_hist_train,'-c 990 -g 0.01 -b 1');

gold_trng_lbl_vector = gold_trng_lbl_vector';
model_gold = svmtrain(gold_trng_lbl_vector,hsv_hist_train,'-c 900 -g 0.01 -b 1');


for k = 1:length(testing_set)  
    testing_label_vector(k) = 0.01;
end
testing_label_vector = testing_label_vector';

[predicted_label_blak, accuracy, prob_estimates_blak] = svmpredict(testing_label_vector, hsv_hist_test, model_blak,'-b 1');
[predicted_label_red, accuracy, prob_estimates_red] = svmpredict(testing_label_vector, hsv_hist_test, model_red,'-b 1');
[predicted_label_silver, accuracy, prob_estimates_silver] = svmpredict(testing_label_vector, hsv_hist_test, model_silver,'-b 1');
[predicted_label_brown, accuracy, prob_estimates_brown] = svmpredict(testing_label_vector, hsv_hist_test, model_brown,'-b 1');
[predicted_label_gold, accuracy, prob_estimates_gold] = svmpredict(testing_label_vector, hsv_hist_test, model_gold,'-b 1');

tempg = prob_estimates_gold(:,2);
[Y,I] = sort(tempg,'descend');
result_gold = {};
for k=1:200
    result_gold{k} = testing_set{I(k)};
end
closestMatches = result_gold;

% create a webpage showing the results
fid = fopen('mining_image_gold.html','w');
fprintf(fid,'<html><body>\n');
for i=1:length(result_gold)
	picname = result_gold{i};
	fprintf(fid,['<img src="bags/' picname  '">\n']);
	if mod(i,3)==0
		fprintf(fid,'<br>');
	end
end
fprintf(fid,'</html>');
fclose(fid);

tempbl = prob_estimates_blak(:,1);
[Y,I] = sort(tempbl,'descend');
result_black = {};
for k=1:200
    result_black{k} = testing_set{I(k)};
end
closestMatches = result_black;

% create a webpage showing the results
fid = fopen('mining_image_black.html','w');
fprintf(fid,'<html><body>\n');
for i=1:length(result_black)
	picname = result_black{i};
	fprintf(fid,['<img src="bags/' picname  '">\n']);
	if mod(i,3)==0  
		fprintf(fid,'<br>');
	end
end
fprintf(fid,'</html>');
fclose(fid);

tempr = prob_estimates_red(:,2);
[Y,I] = sort(tempr,'descend');

result_red = {};
for k=1:200
    result_red{k} = testing_set{I(k)};
end
closestMatches = result_red;

% create a webpage showing the results
fid = fopen('mining_image_red.html','w');
fprintf(fid,'<html><body>\n');
for i=1:length(result_red)
	picname = result_red{i};
	fprintf(fid,['<img src="bags/' picname  '">\n']);
	if mod(i,3)==0
		fprintf(fid,'<br>');
	end
end
fprintf(fid,'</html>');
fclose(fid);

tempbro = prob_estimates_brown(:,2);
[Y,I] = sort(tempbro,'descend');
result_brown = {};
for k=1:200
    result_brown{k} = testing_set{I(k)};
end
closestMatches = result_brown;

% create a webpage showing the results
fid = fopen('mining_image_brown.html','w');
fprintf(fid,'<html><body>\n');
for i=1:length(result_brown)
	picname = result_brown{i};
	fprintf(fid,['<img src="bags/' picname  '">\n']);
	if mod(i,3)==0
		fprintf(fid,'<br>');
	end
end
fprintf(fid,'</html>');
fclose(fid);


tempsilver = prob_estimates_silver(:,2);
[Y,I] = sort(tempsilver,'descend');
result_silver = {};
for k=1:200
    result_silver{k} = testing_set{I(k)};
end


% create a webpage showing the results
fid = fopen('mining_image_silver.html','w');
fprintf(fid,'<html><body>\n');
for i=1:length(result_silver)
	picname = result_silver{i};
	fprintf(fid,['<img src="bags/' picname  '">\n']);
	if mod(i,3)==0
		fprintf(fid,'<br>');
	end
end
fprintf(fid,'</html>');
fclose(fid);


closestMatches = testing_set;




