type
  bucketArray = array of array of Single;

function ComparePixel(p1, p2 : SingleP): Integer;
begin
	if (p1^ = p2^) then begin
		if (Integer(p1) > Integer(p2)) then begin
			Result := 1;
		end else if (Integer(p1) < Integer(p2)) then begin
			Result := -1;
		end else begin
			Result := 0;
		end;
	end else if (p1^ > p2^) then begin
		Result := 1;
	end else if (p1^ < p2^) then begin
		Result := -1;
	end else begin
		Result := 0;
	end;
end;


//Insertion sort for an almost sorted array (arr_index argument means stop+1 so not inclusive)
procedure InsertionSort(var arr: array of single; start, arr_index: integer);
var
  i,j: integer;
  tmp: single;
begin
  for i:=start+1 to arr_index-1 do begin
    if arr[i]- arr[i-1] < 0 then begin //only bother swapping if not ordered properly
			//the following is no longer used as it seems to be slightly slower (since it does a binary search)
      {if (i > 10) then k := BSearch(arr, 0, i, arr[i])
      else begin
        k := i-1;
        //if k >=0 then begin
        while (arr[i] - arr[k] < 0) and (k >= 0) do begin
          k:= k-1;
        end;
        inc(k);
      end;
      j := i;
      while (j > k) do begin
        tmp := arr[j-1];
        arr[j-1] := arr[j];
        arr[j] := tmp;
        j := j-1;
      end;}
			//just swap down until the proper location
      tmp := arr[i];
      j := i;
      while (tmp < arr[j-1]) do begin
        arr[j] := arr[j-1];
        j := j-1;
        if (j = 0) then break;
      end;
      arr[j] := tmp;
    end;
  end;
end;

//Binary Search (index 0 is start and index 1 is end; points are inclusive)
//Important: if the value is not in the array, returns the index of the least value greater than it
function BSearch(var arr: array of Single; index0, index1: Integer;value: Single): Integer;
var
  min, max: Integer;
begin
  Result := -1;
  if index1 < index0 then exit;
  min := index0; max := index1;
  while (min <= max) do begin
    Result := (min+max) div 2;
    if arr[Result]-value > 0 then max := Result -1
    else if arr[Result]-value<0 then min := Result+1
    else break;
  end;
  if arr[Result]<value then inc(Result);
end;

//Computes the correct bucket to place "point" into, this might be able to be optimized...
function GetBucketIndex(point : Single; var bin_limits : array of single; var fit_output : TBoxStats; num_bins : integer): Integer;
var
  z:                Single; //z-score
  index :           Integer; //index of bin to place point into
begin
  //Indexing into the buckets
  z := (point - fit_output.Avg)/fit_output.Sdev;
  if (abs(z) <= 1) then index := trunc(z*83) + (num_bins div 2) //direct index for +/- 1 sigma
  else if z < -1 then index := Bsearch(bin_limits, 0, length(bin_limits) div 6, point) -1 //binary search lower tail
  else index := Bsearch(bin_limits, (length(bin_limits) div 6) *5, length(bin_limits), point)-1; //binary search upper tail
  Result := index;
end;

//Wrapper for pval because indicies start at 0 for median filter routine
function GetPixelValue( CurrentSlot : Integer; i,j : Integer) : Single;
begin
  Result := pval(CurrentSlot,i+1,j+1);
end;

//Removes point from bucket at index
function DeleteFromBucket(point: Single; index : Integer; var arr : bucketArray; var arr_index : array of integer) : Boolean;
var
  i:  integer;
begin
  for i:=0 to arr_index[index]-1 do if arr[index,i] = point then break; //linear search for the correct value to delete
  if i >= arr_index[index] then begin
    Result := False;
  end else begin
    Move(arr[index,i+1],arr[index,i], (arr_index[index]-1-i)*sizeOf(singleP)); //compress bucket
    arr_index[index] := arr_index[index]-1; //decrease size
    Result := True;
  end
end;

//Adds point to bucket at index
function AddToBucket(point: Single; index: Integer; var arr : bucketArray; var arr_index : array of integer) : Boolean;
begin
  arr[index,arr_index[index]] := point;
  inc(arr_index[index]);
end;

// Bucket sort attempt, using a gaussian fit for bucket params
// Should no longer crash for pathological cases
//   Update by Jason Wang -  11/17/2011
// Should now support handling of NaNs and fixed bug when filtering rectangular arrays
// 
function MedianFilterImages(start,stop,dest,box_height,box_width: integer): boolean;
var
  row, col:             integer; //position of center of box
  i,j,k,t1,t2:          integer; //misc variables, time trackers for debugging
  Xmin,Ymin,Xmax,Ymax:  integer; //box size, in absolute indicies
  hlo,hhi,wlo,whi:      integer; //box size, offsets from center
  tmin,tmax:            single; //min, max variables for special case of median
  index, diff:          integer; //variables to index into the array
  arr:                  bucketArray; //the buckets
  arr_index:            array of integer; //array of size of buckets
  num_bins:             integer; //number of buckets
  bin_limits:           array of single; //the limits of each buckets
  box:                  TDrawingTool; //variable for gaussian fitting
  ll, tr:               TPoint; //variables fro gaussian fitting
  fit_params:           TBoundsRecord; //input for gaussian fitting
  fit_output:           TBoxStats; //output for gaussian fitting
  direction:            integer; //direction the box is moving (1 is forward, -1 is backward)
  median:               single; //value for median
  point:                single; //variable to access value of points in the array
  npts:                 integer;
  sum:                  single;
  sn:                   integer;
  ShowProgress:         boolean;
  ShowIntProg:          boolean;
  imin,imax:            integer;
  z:                    single; //z score variable for index calculation
  current_size:         integer; //current size of median filter box
begin
  //t1 := GetTickCount;
  ShowProgress := stop - start + 1 > 2;
  if ShowProgress then DisplayProgress(2,'Median Filter',start,stop);
  //set box size parameters
	wlo := (box_width-1) div 2; hlo := (box_height-1) div 2;
  whi := box_width -1 - wlo; hhi := box_height -1 - hlo;
	//set bucket sizes
  num_bins := 256;
  setlength(bin_limits, num_bins+1); //specifies the edges of the buckets
  setlength(arr_index, num_bins);
  setlength(arr, num_bins);
  for i := 0 to num_bins-1 do setlength(arr[i],box_height*box_width*2);
  for sn := start to stop do begin
    if not Slot[sn].HasData then exit;
    if not AllocateWorkSlot(sn) then exit;
    with Slot[sn] do begin
      //initialize variables
	    row := 0; col := 0;
	    direction := 1;
      for i:=0 to num_bins-1 do arr_index[i]:=0;
      current_size := 0;

      //initialize fit params
      box := dtRectangle;
      ll.X := 1; ll.Y := 1; tr.X:= NCols-1; tr.Y:= NRows-1;
      fit_params.LoReject := false; fit_params.HiReject := false;
      fit_params.TrimLo := true; fit_params.TrimHi := true;
      fit_params.NumLoTrim := 10; fit_params.NumHiTrim := 10;
      fit_params.SDReject := true;
      fit_params.NumSigma := 3.5;
      fit_params.NumIter := 3;
      fit_params.LoLimit := m_infinity;
      fit_params.HiLimit := p_infinity;


      //Compute bucket ranges
      ComputeAreaStats(sn, box, ll, tr, 0, fit_params, fit_output);
      if fit_output.Sdev = 0 then begin Result := true; exit; end; //if uniform image just exit
      bin_limits[num_bins] := 1E38; //nothing should be greater than the max
      bin_limits[0] := -1E38; //nothing should be less than the min
      for i := 1 to length(bin_limits) div 2 do begin
        bin_limits[i] := (-1*z_prob[length(z_prob)-i] * fit_output.Sdev) + fit_output.Avg;
      end;
      for i := length(bin_limits) div 2 + 1 to length(bin_limits) -2 do begin
        bin_limits[i] :=  (z_prob[i - length(bin_limits) div 2] * fit_output.Sdev) + fit_output.Avg;
      end;

      //initialize buckets
      // i = x, j = y
      for i := 0 to whi do begin
	      for j := 0 to hhi do begin
          point := GetPixelValue(sn,i,j);
          if ((point <= m_infinity) OR (point >= p_infinity)) then begin
            continue;
          end else begin
			      //Indexing into the buckets
            index := GetBucketIndex(point, bin_limits, fit_output,num_bins);
            //Put value into bucket, extending the bucket size if necessary
		  			diff := arr_index[index];
            arr[index,diff] := point;
            inc(arr_index[index]);
            current_size := current_size + 1;
          end;
         end;
      end;

      //row = y, col = x
      //begin back and forth pattern
      while row < NRows do begin
        //Set box height parameters
        if (row-hlo >= 0) then Ymin := row-hlo
          else Ymin := 0;
        if (row+hhi < NRows) then Ymax := row+hhi
          else Ymax := NRows-1;

        while ((col < NCols) AND (direction = 1)) OR ((col >=0) AND (direction = -1)) do begin
           //Set box width parameters
          if(col-wlo >= 0) then Xmin := col-wlo
             else Xmin := 0;
          if(col+whi < NCols) then Xmax := col+whi
             else Xmax := NCols-1;

          //Add/delete from left side
          if (col-wlo >= 0) then begin
            //if going forwards, delete from left if needed
            if ((direction > 0) AND (col-wlo <> 0)) then begin
               for j := Ymin to Ymax do begin
                Point := GetPixelValue(sn,Xmin-1,j);
                if ((point <= m_infinity) OR (point >= p_infinity)) then begin
                  continue;
                end else begin
                  index := GetBucketIndex(point, bin_limits, fit_output,num_bins);
                  DeleteFromBucket(point,index,arr,arr_index);
                  current_size := current_size - 1;
                end;
               end;
				    //if going backwards, add if necessary
            end else if ((direction < 0) AND (col <> NCols-1)) then begin
              for j := Ymin to Ymax do begin
                Point := GetPixelValue(sn, Xmin, j);
                if ((point <= m_infinity) OR (point >= p_infinity)) then begin
                  continue;
                end else begin
                  index := GetBucketIndex(point, bin_limits, fit_output,num_bins);
                  AddToBucket(point,index,arr,arr_index);
                  current_size := current_size + 1;
                end;
              end;
            end;
          end;    //  if (col-wlo >= 0) then begin

			    //Add/delete from right side
          if (col+whi < NCols) then begin
            //if going backwards, delete from right if needed
            if ((direction < 0) AND (col+whi <> NCols-1)) then begin
              for j := Ymin to Ymax do begin
                Point := GetPixelValue(sn,Xmax+1,j);
                if ((point <= m_infinity) OR (point >= p_infinity)) then begin
                  continue;
                end else begin
                  index := GetBucketIndex(point, bin_limits, fit_output,num_bins);
                  DeleteFromBucket(point,index,arr,arr_index);
                  current_size := current_size - 1;
                end;
              end;
            //if going forwards, add from right if necessary
            end else if ((direction > 0) AND (col <> 0)) then begin
              for j := Ymin to Ymax do begin
                Point := GetPixelValue(sn, Xmax, j);
                if ((point <= m_infinity) OR (point >= p_infinity)) then begin
                  continue;
                end else begin
                  index := GetBucketIndex(point, bin_limits, fit_output,num_bins);
                  AddToBucket(point,index,arr,arr_index);
                  current_size := current_size + 1;
                end;
              end;
            end;
          end;      //  if (col+whi < NCols) then begin

          //count bin sizes until median
          median := 0;
          index :=0;
          //current_size :=0;
          //for k:=0 to num_bins-1 do current_size := current_size + arr_index[k];//(Ymax-Ymin+1)*(Xmax-Xmin+1)-num_invalid;
          for i :=0 to num_bins-1 do begin
            //if there's nothing in the histogram, keep original value (NaN likely)
            if current_size <= 0 then begin
              median := GetPixelValue(sn, col, row);
              break;
            end;
            index := index + arr_index[i];
            if index > (current_size-1)/2 then begin
              //once count has past the median, look in the last bucket for the median
              if (index = round(current_size/2)) and (current_size mod 2 = 0) then begin
                //unless the count is just on the median and the number of points is even
                //then find max of i bucket and min of the next non-empty bucket after it
                j := i+1;
                while(arr_index[j] < 1) AND (j<num_bins) do j := j+1;
                tmax := Single(arr[i,0]);
                for k:=1 to arr_index[i]-1 do if Single(arr[i,k]) > tmax then tmax := Single(arr[i,k]);
                tmin := Single(arr[j,0]);
                for k:=1 to arr_index[j]-1 do if Single(arr[j,k]) < tmin then tmin := Single(arr[j,k]);
                median := (tmax+tmin)/2;
              end else begin
                //sort array and index into correct location
                InsertionSort(arr[i], 0,arr_index[i]);
                diff := arr_index[i]-(index - (current_size div 2));
                if (current_size MOD 2 = 1) then median := Single(arr[i,diff])
                else median := (Single(arr[i,diff]) + Single(arr[i,diff-1]))/2;
              end;
              break;      //once median is found, end looping
            end;
          end;          //  for i :=0 to num_bins-1 do begin

          set_pval(WorkSlot,col+1,row+1,median); //finally set the median
          col := col + direction; //continue along the row

        end;     //  while ((col < NCols) AND (direction = 1)) OR ((col >=0) AND (direction = -1)) do begin

        direction := -1*direction; //switch directions
        col := col + direction; //at the end, goes off the array, needs to come back in

		    //delete bottom row if necessary
        if (row-hlo >= 0) then begin
          for i := Xmin to Xmax do begin
            Point := GetPixelValue(sn, i, Ymin);
            if ((point <= m_infinity) OR (point >= p_infinity)) then begin
              continue;
            end else begin
              index := GetBucketIndex(point, bin_limits, fit_output,num_bins);
              DeleteFromBucket(point,index,arr,arr_index);
              current_size := current_size - 1;
            end;
          end;
        end;                    //  if (row-hlo >= 0) then begin

		    //add top row if necessary
        if (row+hhi < NRows-1) then begin
          for i := Xmin to Xmax do begin
            Point := GetPixelValue(sn, i, Ymax+1);
            if ((point <= m_infinity) OR (point >= p_infinity)) then begin
              continue;
            end else begin
              index := GetBucketIndex(point, bin_limits, fit_output,num_bins);
              AddToBucket(point,index,arr,arr_index);
              current_size := current_size + 1;
            end;
          end;
        end;                    //  if (row+hhi < NCols-1) then begin

        row := row + 1;                         //go to next line
      end;    //  while row < NRows do begin
    end;     // with Slot[sn] do
    UpdateFITSHdr(Header[WorkSlot],WorkSlot);
    AddFITSHistory(WorkSlot,'Median Filtered: '+IntToStr(box_width)+' x '+
                            IntToStr(box_height));
    Slot[WorkSlot].DispDefined := false;
    CopySlot(WorkSlot,dest);
    dest := dest + 1;
    if ShowProgress and UpdateProgress(2,sn) then begin
      QuitProgress(2,true);
      CloseProgressWindow := true;
      exit;
    end;
  end;    // for sn := start ..
  for i := 0 to length(arr)-1 do arr[i] := nil;
  arr := nil;
  arr_index := nil;
  bin_limits := nil;
  EmptySlot(WorkSlot);
  if ShowProgress then QuitProgress(2,CloseProgressWindow);
  result := true;
  //t2 := GetTickCount;
  //DebugForm.DebugWindow.Output.Lines.Add(format('Time is %5d, time per frame is %6.4gs', [t2-t1, (t2-t1)/(stop-start+1)/1000]));
end;