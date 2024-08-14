function [positions,orientations,objEntries] = readMOT(motionfilename)
    fp = fopen( motionfilename, 'r' );

    tline = fgetl(fp);
    tline = strrep(tline,"NumFrames: ", "");
    numframes = str2num(tline);
    positions = [];
    orientations = [];
    objEntries = {};
   
    while ~feof(fp)
        tline = fgetl(fp);
        if contains(tline,'position')
            bracketNumbers = regexp(tline, '\[(\d+)\]', 'tokens');
            objNameEnd = strfind(tline,'[');
            objName = tline(1:objNameEnd-1);
            objIndex = str2num(bracketNumbers{1}{1});
            lastEntry = numel(objEntries);
            objEntries{lastEntry+1} = {objName,objIndex};

            objPositions = zeros(numframes,3);
            for i = 1:numframes
                tline = fgetl(fp);
                nums = str2num(tline);
                objPositions(i,:) = nums;
            end
            positions = cat(3,positions,objPositions);
            continue;
        end

        if contains(tline,'orientation') %assumes the orientation is right after the position
            bracketNumbers = regexp(tline, '\[(\d+)\]', 'tokens');
            objIndex = str2num(bracketNumbers{1}{1});

            objOrientations = zeros(numframes,4);
            for i = 1:numframes
                tline = fgetl(fp);
                nums = str2num(tline);
                objOrientations(i,:) = nums;
            end
            orientations = cat(3,orientations,objOrientations);
            continue;
        end
    end

    fclose(fp);
end

