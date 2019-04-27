#include "DP_matrix.h"

using namespace std;

/*
Computes the mutation costs according to the mut-function of the theoretical framework
*/
float mut(int refA, int refB, int hA, int hB, float c) {
	return c*(abs((refA-hA))+abs((refB-hB)));
}

/*
	performs the actual computation of the scoring matrix and returns this
*/
tuple<float**, vector<int>> perform_DP(vector<int> H_A, vector<int> H_B, vector<vector<int>> Ref, float param_mut, unordered_set<int> E){
	int m = H_A.size();
	int n = Ref.size();
	cout << "m: " << m << endl;
	cout << "n: " << n << endl;
	//Initialize scoring matrix	
	int cols = ceil(sqrt(m));
//	int cols = floor(sqrt(m));
	float** Score = new float*[cols];
	for(int i = 0; i < cols+1; ++i)
		Score[i] = new float[n*n];
	cout << "number of columns: " << cols << endl;	
	//Initialize the arrays containing minima of columns
//	float* minrightlist = new float[n];
//	float* minleftlist = new float[n];
	vector<float> minrightlist;
	vector<float> minleftlist;
	for (int k = 0; k < n; ++k) {
		minrightlist.push_back(INT_MAX);		
		minleftlist.push_back(INT_MAX);	
	}
	cout << "length minrightlist: "<< minrightlist.size() << endl;
	cout << "minrightlist 0: "<< minrightlist[0] << endl;
	
	int x = Ref.size();
	int y = Ref[0].size();
	
	//initialize the reference matrix
	int** ref = new int*[x];
	for(int i = 0; i < x; ++i)
		ref[i] = new int[y];
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
		 	ref[i][j] = Ref[i][j];	
		}
	}
	
	//converts the haplotype vectors to int arrays
	int* ha = new int[H_A.size()];
	int* hb = new int[H_B.size()];
	for (unsigned int i = 0; i < H_A.size(); i++) {
		ha[i] = H_A[i];	
		hb[i] = H_B[i];
	}

	float mingeneral = INT_MAX;
	//initialize first row of scoring matrix
	for(int left = 0; left < n; left++) {
		for (int right = 0; right < n; right++)  {
			int j = left*n + right;
			Score[0][j] = mut(Ref[left][0],Ref[right][0],H_A[0],H_B[0],param_mut);
			if (Score[0][left*n + right] < minrightlist[right]){
				minrightlist[right] = Score[0][left*n + right];
			}
			if (Score[0][left*n + right] < minleftlist[left]){
				minleftlist[left] = Score[0][left*n + right];
			}
			if (Score[0][j] < mingeneral) mingeneral = Score[0][left*n + right];
		}
	}
	cout << "Score 0 0 before: " << Score[0][0] << endl;
	vector<int> stored_positions;
	stored_positions.push_back(0);
	float  mutcost, switchcost, score = 0;
	float c00, c01, c10,c11,c02,c12,c21 = 0;
	int j = 0;
	vector<float> temp_column;
	for (int i = 1; i < m; ++i) {
		vector<float> newminrightlist;
		vector<float> newminleftlist;
		for (int k = 0; k < n; ++k) {
			newminrightlist.push_back(INT_MAX);	
			newminleftlist.push_back(INT_MAX);	
		}		
		cout << "i: " << i << endl;
	//	if (i%cols == 0 || i==m-1) 
		if ((i+1) % (cols-1) == 0 || i == m-1)
			stored_positions.push_back(i);
		//cout << "i%cols: " << i%cols << endl;
		//fills the scoring matrix using the efficient computation of  O//(m^2 x n)
		for (int left = 0; left < n; left++) {
		//	cout << "left: " << left << endl;
			for (int right = 0; right < n; right++)  {
			//	cout << "rigth: " << right << endl;
				j = left*n + right;
			//	cout << "j: " << j << endl;
				mutcost = mut(ref[left][i],ref[right][i],ha[i],hb[i], param_mut);
			//	cout << "mutcost: " << mutcost << endl;
				if (i == 1){
					c00 = Score[i-1][j];
				}
				else {
					//todo: is this index j correct here?
			//		cout << "temp column in else: " << temp_column[j] << endl;
					c00 = temp_column[j];
				}
				c01 = 1+minrightlist[right];
			//	cout << "c01 "<< c01 << endl;				
				c10 = 1+minleftlist[left];
				c11 = 2+mingeneral;
				
				if (E.find(i-1) != E.end()) {
					//todo: is this index correct here?
					c02 = temp_column[right*n + left];
					c12 = 1+minrightlist[left];
					c21 = 1+minleftlist[right];
					float mins[] = {c00,c01,c10,c11,c02,c12,c21};			
					switchcost = *min_element(mins,mins+7);
				}
				else {
					float mins[] = {c00,c01,c10,c11};
					switchcost = *min_element(mins, mins+4);				
				}
				score = switchcost + mutcost;
				if ((i+1) % (cols-1) == 0 || i == m-1) {
					//stored_positions.push_back(i);
					//Score[(int)(i/cols)][j] = switchcost + mutcost;
					Score[(int)((i+1)/(cols-1))][j] = switchcost + mutcost;
				}
			//	cout << "Score: " << score << endl;
				if (i == 1) temp_column.push_back(score);
				else temp_column.at(j) = score;				
				if (score < newminrightlist[right]) newminrightlist[right] = score;				
				if (score < newminleftlist[left]) newminleftlist[left] = score;				
				if (score < mingeneral) mingeneral = score;		
			}
		}
		mingeneral = INT_MAX;
		minrightlist = newminrightlist;
		minleftlist = newminleftlist;
	}
//	delete [] minrightlist;
//	delete [] minleftlist;
	cout << "Score[0]0 " << Score[0][0] << endl;
	cout << "Score[1]0 " << Score[1][0] << endl;
	
	for(int i = 0; i < x; ++i) {
		delete [] ref[i];
	}
	delete [] ref;
	delete [] ha;
	delete [] hb;
	return make_tuple(Score, stored_positions );
	for(int i = 0; i < cols+1; ++i) {
    delete [] Score[i];
	}
	delete [] Score;
	
	
}	

tuple<float**,int> DP_partial(vector<int> H_A, vector<int> H_B, vector<vector<int>> Ref, float param_mut, int start, int end, float** complete_Score, unordered_set<int> E, map<int,int> complete_to_partial_index) {
	//compute DP between Score[start] and Score[end]
	int n = Ref.size();

	//Initialize scoring matrix	
	int cols = end-start+1;
	cout << "in DP_partial: start: " << start << endl;
	cout << "in DP_partial: end: " << end << endl;
	cout << "in DP_partial: cols: " << cols << endl;
	float** Score = new float*[cols];
	for(int i = 0; i < cols; ++i)
		Score[i] = new float[n*n];

	vector<float> minrightlist;
	vector<float> minleftlist;
	for (int k = 0; k < n; ++k) {
		minrightlist.push_back(INT_MAX);	
		minleftlist.push_back(INT_MAX);	
	}
	
	int x = Ref.size();
	int y = Ref[0].size();
	
	//initialize the reference matrix
	int** ref = new int*[x];
	for (int i = 0; i < x; ++i)
		ref[i] = new int[y];
	for (int i = 0; i < x; i++) {
		for (int j = 0; j < y; j++) {
		 	ref[i][j] = Ref[i][j];	
		}
	}
	
	//converts the haplotype vectors to int arrays
	int* ha = new int[H_A.size()];
	int* hb = new int[H_B.size()];
	for (unsigned int i = 0; i < H_A.size(); i++) {
		ha[i] = H_A[i];	
		hb[i] = H_B[i];
	}
	
	int newstart = complete_to_partial_index[start];
	int newend = complete_to_partial_index[end];	
	float* startrow = complete_Score[newstart];
	float* endrow = complete_Score[newend];
	cout << "in DP partial, newstart: " << newstart << endl;
	cout << "in DP partial, newend: " << newend << endl;
	cout << "in DP partial, startrow[0]: " << startrow[0] << endl;
	float mingeneral = INT_MAX;
	//initialize first row of scoring matrix
	for(int left = 0; left < n; left++) {
		for (int right = 0; right < n; right++)  {
			int j = left*n + right;
			Score[0][j] = startrow[j];
			if (startrow[j] < minrightlist[right]){
				minrightlist[right] = startrow[j];
			}
			if (startrow[j] < minleftlist[left]){
				minleftlist[left] = startrow[j];
			}		
			if (startrow[j] < mingeneral) mingeneral = startrow[j];
		}
	}
	
	float  mutcost, switchcost, score = 0;
	float c00, c01, c10,c11,c02,c12,c21 = 0;
	int j = 0;

	for (int i = start+1; i < end+1; ++i) {
		cout << "in DP partial, i: " << i << endl;
		vector<float> newminrightlist;
		vector<float> newminleftlist;
		for (int k = 0; k < n; ++k) {
			newminrightlist.push_back(INT_MAX);	
			newminleftlist.push_back(INT_MAX);	
		}
		//fills the scoring matrix using the efficient computation of  O(m^2 x n)
		for (int left = 0; left < n; left++) {
			for (int right = 0; right < n; right++)  {
				j = left*n + right;
					
				mutcost = mut(ref[left][i],ref[right][i],ha[i],hb[i], param_mut);
				c00 = Score[i-start-1][j];
				c01 = 1+minrightlist[right];
				c10 = 1+minleftlist[left];
				c11 = 2+mingeneral;
				if (E.find(i-1) != E.end()) {
					//todo: is this index correct here?
					c02 = Score[i-start-1][right*n + left];
					c12 = 1+minrightlist[left];
					c21 = 1+minleftlist[right];
					float mins[] = {c00,c01,c10,c11,c02,c12,c21};			
					switchcost = *min_element(mins,mins+7);
				}
				else {
					float mins[] = {c00,c01,c10,c11};
					switchcost = *min_element(mins, mins+4);				
				}
				score = switchcost + mutcost;
				Score[i-start][j] = score;
				
				if (score < newminrightlist[right]) newminrightlist[right] = score;				
				if (score < newminleftlist[left]) newminleftlist[left] = score;				
				if (score < mingeneral) mingeneral = score;

			//	}			
			}
		}
		mingeneral = INT_MAX;
		minrightlist = newminrightlist;
		minleftlist = newminleftlist;
	}
	cout << "in DP partial, cols: " << cols << endl;
	
	for (int i = 0; i < x; ++i)
		delete [] ref[i];
	delete [] ref;
	delete [] ha;
	delete [] hb;
	return make_tuple(Score, cols);
	for(int i = 0; i < cols; ++i) {
    delete [] Score[i];
	}
	delete [] Score;
}

tuple<float, int> minimum(float* current_list, int l, int r, bool is_terminal, int N) {
	/*
	Needed for the backtracing. Performs for a given column in current_list the backward computation
	 of the scoring function to determine the optimal value it origins from. 
	 Output: The minimum value in the preceding column and the according index
	*/
	vector<float> new_list;
	int size = N*N;
	int n = (int)sqrt(size);
	int switchval = 0;
	//backtracing from one column in the DP matrix to a column that was at the end of a block
	if (is_terminal) {
		for (int i = 0; i < size; i++) {			
			int l_ = (int)(i/n);
			int r_ = (int)(i%n);
			if ((l_ == r && r_==l) || (l_ == l && r_ == r))
				switchval = 0;
			else if ((l_ == r && r_ != l) || (l_ != r && r_ == l))
				switchval = 1;
			else if ((l_ == l && r_ != r) || (l_ != l && r_ == r))
				switchval = 1;
			else
				switchval = 2;
			new_list.push_back(current_list[i]+switchval);
		}
		//find optimal value when considering switch costs between columns
		float newmin = *min_element(new_list.begin(), new_list.end());
		int newmin_index = find(new_list.begin(), new_list.end(), newmin) - new_list.begin();
		return make_tuple(newmin, newmin_index);
	}
	//backtracing to a column that was not the end of a block
	else {
		for (int i = 0; i < size; i++) {
			int l_ = (int)(i/n);
			int r_ = (int)(i%n);
			if (l_ != l && r_ !=r)
				switchval = 2;
			else if (r_ != r && l_ == l)
				switchval = 1;
			else if (r_ == r && l_ != l)
				switchval = 1;
			else if (r_ == r && l_ == l) switchval = 0;
			new_list.push_back(current_list[i]+switchval);
		}
		int newmin = *min_element(new_list.begin(), new_list.end());
		int newmin_index = find(new_list.begin(), new_list.end(), newmin) - new_list.begin();
		return make_tuple(newmin, newmin_index);
	}		
}

tuple<vector<int>,vector<int>> backtrace_partial(float** Score, int endpos, int n, float** Score_partial, unordered_set<int> e, map<int,int> complete_to_partial_index, int cols) {
	//determines the starting point for the backtracing
	vector<int> path1;
	vector<int> path2;

	cout  << "in backtrace, endpos: " << endpos << endl;
	cout << "in backtrace, complete_to_partial_index[endpos]: " << complete_to_partial_index[endpos] << endl;
	float mini = Score[complete_to_partial_index[endpos]][0];
	cout << "in backtrace:, mini " << mini << endl;
	int minarg = 0;
	int l = 0;
	int r = 0;
	for (int i = 0; i < n*n; i++) {
		if (Score[complete_to_partial_index[endpos]][i] < mini) {
			minarg = i;		
			mini = Score[complete_to_partial_index[endpos]][i];
		}
	}
	cout << "in backtrace_partial, minarg: " << minarg << endl;
	//loop over Score_partial
//	int end = Score_partial.size() -1;
	int end = (sizeof(Score_partial)/sizeof(*Score_partial))-1;
	//end = cols;
	//as a test
/*	for (int i = 0; i < n*n; i++) {
		if (Score_partial[end][i] < mini) {
			minarg = i;		
			mini = Score_partial[end][i];
		}
	}
*/
	cout << "in backtrace_partial, minarg 2: (should be equal minarg before) " << minarg << endl;
	l = (int)minarg/n;
	r = (int)minarg%n;
	path1.push_back(l);
	path2.push_back(r);

	cout << "complete_to_partial_index[endpos] (should equal end?) " << complete_to_partial_index[endpos] << endl;
	cout << "end: " << end << endl;
	end = cols;
	cout << "end after: " << end << endl;
	while(end > 0) {
		minarg = get<1>(minimum(Score_partial[end-1], l, r, e.find(endpos)!= e.end(), n) );
		l = (int)minarg/n;
		r = (int)minarg%n;

		path1.push_back(l);
		path2.push_back(r);
		end -= 1;
		endpos-= 1;
		//todo: test whether complete_to_partial_index[endpos] == end holds in every iteration
	}
	reverse(path1.begin(),path1.end());
	reverse(path2.begin(),path2.end());
	cout << "path1 size: " << path1.size() << endl;	
	return make_tuple(path1, path2);
}

/*
	Reads the input and calls the function for computation of the scoring matrix
*/
void compute_scoring(char* haplofile, char* E_file, char* panelfile, float mutation, char* out_costfile, char* out_pathfile){	
	//haplofile: A file with the haplotype strings extracted from the target vcf
	
	ifstream file(haplofile);
 	string haplo1;
 	string haplo2;
	getline(file, haplo1);
   getline(file, haplo2);
   
   //E_file: A file containing the indices for the block ending boundaries
  /* int counter =-1;
 		
 	ifstream infile(E_file);
  	string line;
	while(getline(infile,line)){
   	stringstream linestream(line);
    	string value;
   	while(getline(linestream,value,',')){
   		counter += 1;
		}
	}
	*/
	vector<int> vectorE;
	ifstream infile2(E_file);
	int counter2 = -1;
  	string line2;
	while(getline(infile2,line2)){
   	stringstream linestream(line2);
    	string value;
   	while(getline(linestream,value,',')){
   		counter2 += 1;
      	vectorE.push_back(stoi(value));
		}
	}
	unordered_set<int> e{begin(vectorE), end(vectorE)};
	
	//push strings for the haplotypes into vectors of ints
	vector<int> H_A;
	for (unsigned int i = 0; i < haplo1.size(); i++) {
		H_A.push_back(haplo1[i] - '0');	
	}
	vector<int> H_B;
	for (unsigned int i = 0; i < haplo2.size(); i++) {
		H_B.push_back(haplo2[i] - '0');	
	}
	
	//panelfile: The txt file containing the formerly created reference panel
	ifstream ref_panel(panelfile);
	string newline;
	//create a reference matrix
	vector<vector<int> > Ref;
	while (getline(ref_panel, newline)) {
		istringstream iss(newline);
		char ch;
		vector<int> row;
		while(iss >> ch) {
			row.push_back(ch-'0');		
		}
		Ref.push_back(row);
	}

	//read in the mutation parameter
	float param_mut = mutation;
	
	//perform the actual computation of the scoring matrix and return the matrix 'Score'
	float** Score;
	vector<int> stored_positions;	
	tie(Score,stored_positions) = perform_DP(H_A,H_B,Ref, param_mut,e);
	cout << "stored positions size: " << stored_positions.size() <<  endl;
	cout << "stored positions at 0: " << stored_positions.at(0) <<  endl;
	cout << "stored positions at 1: " << stored_positions.at(1) <<  endl;
	cout << "stored positions at 174: " << stored_positions.at(174) <<  endl;
	cout << "stored positions back: " << stored_positions.back() <<  endl;
	map<int,int> complete_to_partial_index;
	for (int i = 0; i < stored_positions.size(); i++) {
		complete_to_partial_index[stored_positions[i]] = i;	
	}
	cout << "Score[0]: " << Score[0] << endl;
	
	int end = H_A.size() -1;
	cout << "end: " << end << endl;
	int last_pos = stored_positions.at(stored_positions.size()-1);
	//todo: check whether these two are the same
	cout << "last pos: (should equal end) " << last_pos << endl;
	
	int n = Ref.size();

	//determines the starting point for the backtracing
	vector<int> path1;
	vector<int> path2;
	float mini = Score[end][0];
	int minarg = 0;
	int l = 0;
	int r = 0;
	for (int i = 0; i < n*n; i++) {
		if (Score[end][i] < mini) {
			minarg = i;		
			mini = Score[end][i];
		}
	}

	//backtracing to create the paths through the scoring matrix
	l = (int)minarg/n;
	r = (int)minarg%n;
//	path1.push_back(l);
//	path2.push_back(r);
	int stored_size = stored_positions.size();
	float** Score_partial;	
	int cols;
	vector<int> part1;
	vector<int> part2;
	while (stored_size > 1) {
	/*	minarg = get<1>(minimum(Score[end-1], l, r, e.find(end-1)!= e.end(), n) );
		l = (int)minarg/n;
		r = (int)minarg%n;

		path1.push_back(l);
		path2.push_back(r);
		end -= 1;
	*/	
		int startpos = stored_positions.at(stored_size-2);
		int endpos = stored_positions.at(stored_size-1);
		cout << "startpos: " << startpos << endl;
		cout << "endpos: " << endpos << endl;
	//	cout << "pos before startpos: " << stored_positions.at(stored_size-3) << endl;
		// recompute part between Score[stored_positions.at(stored_positions.size()-2)] and Score[stored_positions.at(stored_positions.size()-1)]: Score = dp_partial(H_A, H_B, Ref, param_mut, end-2, end-1, Score, e, complete_to_partial_index)
			
		tie(Score_partial, cols) = DP_partial(H_A, H_B, Ref, param_mut, startpos, endpos, Score, e, complete_to_partial_index);
		// perform backtracing in this part, beginning: l and r are given through already computed minarg
		part1.clear();
		part2.clear();
		cout <<" in compute scoring, cols: " << cols << endl;
		cout << "DP partial finished, now backtrace partial " << endl;
		cout << "Score[174][0]: " << Score[174][0] << endl;
		cout << "Score[174][100]: " << Score[174][100] << endl;
		cout << "Score[173][0]: " << Score[173][0] << endl;
		cout << "Score[173][100]: " << Score[173][100] << endl;
		cout << "Score[0][100]: " << Score[0][100] << endl;	
		cout << "Score[175][0]: " << Score[175][0] << endl;
		cout << "Score[176][0]: " << Score[176][0] << endl;
		tie(part1, part2) = backtrace_partial(Score, endpos, n, Score_partial, e, complete_to_partial_index, cols);
		//  assembled path parts
		path1.insert(path1.begin(), part1.begin()+1, part1.end()-1);
		path2.insert(path2.begin(), part2.begin()+1, part2.end()-1);
	cout << "path1 length: " << path1.size() << endl;
		stored_size -= 1;
		for(int i = 0; i < cols; ++i) {
			delete [] Score_partial[i];
		}
		delete [] Score_partial;
	}
//	for (int i =0; i < path1.size(); i++) {
//		cout << "path: " << path1[i] << endl;	
//	}

	//paths were assembled from back, so they need to be reversed
//	reverse(path1.begin(),path1.end());
//	reverse(path2.begin(),path2.end());

	//outfile for the paths 	
	ofstream myfile;
   myfile.open(out_pathfile);
	for (auto &i : path1) {
		myfile << i << " ";
	}
	myfile << '\n';
	for (auto &i : path2) {
		myfile << i << ' ';
	}
	myfile << '\n';
}