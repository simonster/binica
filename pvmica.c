/******************************************************************************/
/* For all support and information please contact:                            */
/*                                                                            */
/*   Sigurd Enghoff                                                           */
/*   The Salk Institute, CNL                                                  */
/*   enghoff@salk.edu                                                         */
/*                                                                            */
/* Additional ICA software:                                                   */
/*   http://www.cnl.salk.edu/~enghoff/                                        */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "pvmica.h"


/********************* Concatenate binary float matricies *********************/
/* Convert size doublereal elements of mat to floating point values and       */
/* appends values to the binary file specified by file.                       */
/*                                                                            */
/* file: FILE pointer (input)                                                 */
/* size: int (input)                                                          */
/* mat:  doublereal array [size] (input)                                      */

void fbc_matwrite(FILE *file, int size, doublereal *mat) {
	float *buffer;
	int i, items;

	if (!file) error("open failed");

	buffer = (float*)malloc(size*sizeof(float));

	for (i=0 ; i<size ; i++) buffer[i] = (float)mat[i];
	
	items = (int)fwrite(buffer,sizeof(float),size,file);
	if (items != size) error("invalid number of elements");
	
	free(buffer);
}


/********************* Concatenate ascii integer matricies ********************/
/* Append size integer values of mat as ascii numbers to the file specified   */
/* by file.                                                                   */
/*                                                                            */
/* file: FILE pointer (input)                                                 */
/* size: int (input)                                                          */
/* mat:  integer array [size] (input)                                         */

void iac_matwrite(FILE *file, int size, integer *mat) {
	int i;

	if (!file) error("open failed");
	for (i=0 ; i<size-1 ; i++) fprintf(file,"%d ",(int)(mat[i]));
	fprintf(file,"%d\n",(int)(mat[i]));
}


/************************** Send assignment to slave **************************/
/* Send assignment structure assign to slave process tid.                     */
/*                                                                            */
/* tid:    int (input)                                                        */
/* assign: tassign pointer (input)                                            */

void send_assign(int tid, tassign *assign) {
	int     bufid, datasize, i, toint[ASSIGN_NINT];
	double *tdbl, todbl[ASSIGN_NDBL];
	float  *tflt;

	toint[ASSIGN_ID]         = (int)(assign->id);
	toint[ASSIGN_CHANS]      = (int)(assign->chans);
	toint[ASSIGN_FRAMES]     = (int)(assign->frames);
	toint[ASSIGN_EPOCHS]     = (int)(assign->epochs);
	toint[ASSIGN_BIAS]       = (int)(assign->bias);
	toint[ASSIGN_SIGNS]      = (int)(assign->signs);
	toint[ASSIGN_EXTENDED]   = (int)(assign->extended);
	toint[ASSIGN_EXTBLOCKS]  = (int)(assign->extblocks);
	toint[ASSIGN_PDFSIZE]    = (int)(assign->pdfsize);
	toint[ASSIGN_NSUB]       = (int)(assign->nsub);
	toint[ASSIGN_VERBOSE]    = (int)(assign->verbose);
	toint[ASSIGN_BLOCK]      = (int)(assign->block);
	toint[ASSIGN_MAXSTEPS]   = (int)(assign->maxsteps);
	
	todbl[ASSIGN_LRATE]      = (double)(assign->lrate);
	todbl[ASSIGN_ANNEALSTEP] = (double)(assign->annealstep);
	todbl[ASSIGN_ANNEALDEG]  = (double)(assign->annealdeg);
	todbl[ASSIGN_NOCHANGE]   = (double)(assign->nochange);
	todbl[ASSIGN_MOMENTUM]   = (double)(assign->momentum);

	bufid = pvm_initsend(PvmDataDefault);
	pvm_pkint(toint,ASSIGN_NINT,1);
	pvm_pkdouble(todbl,ASSIGN_NDBL,1);

	datasize = assign->chans * assign->frames * assign->epochs;
	if (sizeof(float) != sizeof(doublereal)) {
		tflt = (float*)malloc(datasize*sizeof(float));
		for (i=0 ; i<datasize ; i++) tflt[i] = (float)(assign->data[i]);
		pvm_pkfloat(tflt,datasize,1);
		free(tflt);
	}
	else
		pvm_pkfloat((float*)(assign->data),datasize,1);
		
	datasize = assign->chans * assign->chans;
	if (sizeof(double) != sizeof(doublereal)) {
		tdbl = (double*)malloc(datasize*sizeof(double));
		for (i=0 ; i<datasize ; i++) tdbl[i] = (double)(assign->weights[i]);
		pvm_pkdouble(tdbl,datasize,1);
		free(tdbl);
	}
	else
		pvm_pkdouble((double*)(assign->weights),datasize,1);

	pvm_send(tid,1);
	pvm_freebuf(bufid);
}


/********************************* Kill slave *********************************/
/* Send kill signal to slave process tid.                                     */
/*                                                                            */
/* tid: int (input)                                                           */

void send_kill(int tid) {
	int bufid;

	bufid = pvm_initsend(PvmDataDefault);
	pvm_send(tid,2);
	pvm_freebuf(bufid);
}


/*********************** Receive assignment from master ***********************/
/* Receive assignment structure assign from master.                           */
/*                                                                            */
/* assign: tassign pointer (output)                                           */

int receive_assign(tassign *assign) {
	int     bufid, tag, dummy, datasize, i, toint[ASSIGN_NINT];
	double *tdbl, todbl[ASSIGN_NDBL];
	float  *tflt;

	bufid = pvm_recv(-1,-1);
	pvm_bufinfo(bufid,&dummy,&tag,&dummy);
	if (tag == 2) return 0;
	
	pvm_upkint(toint,ASSIGN_NINT,1);
	pvm_upkdouble(todbl,ASSIGN_NDBL,1);

	assign->id         = (integer)toint[ASSIGN_ID];
	assign->chans      = (integer)toint[ASSIGN_CHANS];
	assign->frames     = (integer)toint[ASSIGN_FRAMES];
	assign->epochs     = (integer)toint[ASSIGN_EPOCHS];
	assign->bias       = (integer)toint[ASSIGN_BIAS];
	assign->signs      = (integer)toint[ASSIGN_SIGNS];
	assign->extended   = (integer)toint[ASSIGN_EXTENDED];
	assign->extblocks  = (integer)toint[ASSIGN_EXTBLOCKS];
	assign->pdfsize    = (integer)toint[ASSIGN_PDFSIZE];
	assign->nsub       = (integer)toint[ASSIGN_NSUB];
	assign->verbose    = (integer)toint[ASSIGN_VERBOSE];
	assign->block      = (integer)toint[ASSIGN_BLOCK];
	assign->maxsteps   = (integer)toint[ASSIGN_MAXSTEPS];
	
	assign->lrate      = (doublereal)todbl[ASSIGN_LRATE];
	assign->annealstep = (doublereal)todbl[ASSIGN_ANNEALSTEP];
	assign->annealdeg  = (doublereal)todbl[ASSIGN_ANNEALDEG];
	assign->nochange   = (doublereal)todbl[ASSIGN_NOCHANGE];
	assign->momentum   = (doublereal)todbl[ASSIGN_MOMENTUM];

	datasize = assign->chans * assign->frames * assign->epochs;
	assign->data = (doublereal*)malloc(datasize*sizeof(doublereal));
	if (sizeof(float) != sizeof(doublereal)) {
		tflt = (float*)malloc(datasize*sizeof(float));
		pvm_upkfloat(tflt,datasize,1);
		for (i=0 ; i<datasize ; i++) assign->data[i] = (doublereal)tflt[i];
		free(tflt);
	}
	else
		pvm_upkfloat((float*)(assign->data),datasize,1);

	datasize = assign->chans * assign->chans;
	assign->weights = (doublereal*)malloc(datasize*sizeof(doublereal));
	if (sizeof(double) != sizeof(doublereal)) {
		tdbl = (double*)malloc(datasize*sizeof(double));
		pvm_upkdouble(tdbl,datasize,1);
		for (i=0 ; i<datasize ; i++) assign->weights[i] = (doublereal)tdbl[i];
		free(tdbl);
	}
	else
		pvm_upkdouble((double*)(assign->weights),datasize,1);

	pvm_freebuf(bufid);
	return 1;
}


/*************************** Send results to master ***************************/
/* Send result structure result to master process tid.                        */
/*                                                                            */
/* tid:    int (input)                                                        */
/* result: tresult pointer (input)                                            */

void send_result(int tid, tresult *result) {
	int    *tint, bufid, datasize, i, toint[RESULT_NINT];
	double *tdbl, todbl[RESULT_NDBL];

	toint[RESULT_ID]    = (int)(result->id);
	toint[RESULT_CHANS] = (int)(result->chans);
	toint[RESULT_BIAS]  = (int)(result->bias != NULL);
	toint[RESULT_SIGNS] = (int)(result->signs != NULL);

	todbl[RESULT_LRATE] = (double)(result->lrate);
	
	bufid = pvm_initsend(PvmDataDefault);
	pvm_pkint(toint,RESULT_NINT,1);
	pvm_pkdouble(todbl,RESULT_NDBL,1);

	datasize = result->chans * result->chans;
	if (sizeof(double) != sizeof(doublereal)) {
		tdbl = (double*)malloc(datasize*sizeof(double));
		for (i=0 ; i<datasize ; i++) tdbl[i] = (double)(result->weights[i]);
		pvm_pkdouble(tdbl,datasize,1);
		free(tdbl);
	}
	else
		pvm_pkdouble((double*)(result->weights),datasize,1);

	if (result->bias != NULL) {
		datasize = result->chans;
		if (sizeof(double) != sizeof(doublereal)) {
			tdbl = (double*)malloc(datasize*sizeof(double));
			for (i=0 ; i<datasize ; i++) tdbl[i] = (double)(result->bias[i]);
			pvm_pkdouble(tdbl,datasize,1);
			free(tdbl);
		}
		else
			pvm_pkdouble((double*)(result->bias),datasize,1);
	}

	if (result->signs != NULL) {
		datasize = result->chans;
		if (sizeof(int) != sizeof(integer)) {
			tint = (int*)malloc(datasize*sizeof(int));
			for (i=0 ; i<datasize ; i++) tint[i] = (int)(result->signs[i]);
			pvm_pkint(tint,datasize,1);
			free(tint);
		}
		else
			pvm_pkint((int*)(result->signs),datasize,1);
	}
	
	pvm_send(tid,1);
	pvm_freebuf(bufid);
}


/************************** Receive result from slave *************************/
/* Receive result structure result from the first slave to responding.        */
/*                                                                            */
/* result: tresult pointer (output)                                           */

int receive_result(tresult *result) {
	int    *tint, bufid, datasize, bias, signs, toint[RESULT_NINT];
	int     i, buflen, tag, tid;
	double *tdbl, todbl[RESULT_NDBL];

	bufid = pvm_recv(-1,-1);
	pvm_upkint(toint,RESULT_NINT,1);
	pvm_upkdouble(todbl,RESULT_NDBL,1);

	result->id         = (integer)toint[RESULT_ID];
	result->chans      = (integer)toint[RESULT_CHANS];
	bias               = (integer)toint[RESULT_BIAS];
	signs              = (integer)toint[RESULT_SIGNS];

	result->lrate      = (doublereal)todbl[RESULT_LRATE];

	datasize = result->chans * result->chans;
	result->weights = (doublereal*)malloc(datasize*sizeof(doublereal));
	if (sizeof(double) != sizeof(doublereal)) {
		tdbl = (double*)malloc(datasize*sizeof(double));
		pvm_upkdouble(tdbl,datasize,1);
		for (i=0 ; i<datasize ; i++) result->weights[i] = (doublereal)tdbl[i];
		free(tdbl);
	}
	else
		pvm_upkdouble((double*)(result->weights),datasize,1);

	if (bias) {
		datasize = result->chans;
		result->bias = (doublereal*)malloc(datasize*sizeof(doublereal));
		if (sizeof(double) != sizeof(doublereal)) {
			tdbl = (double*)malloc(datasize*sizeof(double));
			pvm_upkdouble(tdbl,datasize,1);
			for (i=0 ; i<datasize ; i++) result->bias[i] = (doublereal)tdbl[i];
			free(tdbl);
		}
		else
			pvm_upkdouble((double*)(result->bias),datasize,1);
	}
	else
		result->bias = NULL;

	if (signs) {
		datasize = result->chans;
		result->signs = (integer*)malloc(datasize*sizeof(integer));
		if (sizeof(int) != sizeof(integer)) {
			tint = (int*)malloc(datasize*sizeof(int));
			pvm_upkint(tint,datasize,1);
			for (i=0 ; i<datasize ; i++) result->signs[i] = (integer)tint[i];
			free(tint);
		}
		else
			pvm_upkint((int*)(result->signs),datasize,1);
	}
	else
		result->signs = NULL;

	pvm_bufinfo(bufid,&buflen,&tag,&tid);
	pvm_freebuf(bufid);
	
	return tid;
}


/***************************** Main slave routine *****************************/
/* Repeatedly performs ICA on data segments. Assignments are received by the  */
/* process which decomposes transmitted data and returns the results to the   */
/* assigning master process. The sequence is repeated until a kill signal is  */
/* received.                                                                  */

int slave() {
	int         i, datasize, id, ptid;
	tassign     assign;
	tresult     result;
	doublereal *bias, *data, *weights;
	integer    *signs, chans, frames, epochs;

#ifndef CRAY
	setpriority(PRIO_PROCESS,0,SLAVE_NICE);
#endif

	while (receive_assign(&assign)) {
		id         = assign.id;
		data       = assign.data;
		weights    = assign.weights;
		chans      = assign.chans;
		frames     = assign.frames;
		epochs     = assign.epochs;

		extended   = assign.extended;
		extblocks  = assign.extblocks;
		pdfsize    = assign.pdfsize;
		nsub       = assign.nsub;
		verbose    = assign.verbose;
		block      = assign.block;
		maxsteps   = assign.maxsteps;
		lrate      = assign.lrate;
		annealstep = assign.annealstep;
		annealdeg  = assign.annealdeg;
		nochange   = assign.nochange;
		momentum   = assign.momentum;

		if (assign.bias)
			bias = (doublereal*)malloc(chans*sizeof(doublereal));
		else
			bias = NULL;

		if (assign.signs)
			signs = (integer*)malloc(chans*sizeof(integer));
		else
			signs = NULL;

		runica(data,weights,chans,frames,epochs,bias,signs);
	
		result.id      = id;
		result.weights = weights;
		result.chans   = chans;
		result.bias    = bias;
		result.signs   = signs;
		result.lrate   = lrate;

		ptid = pvm_parent();
		send_result(ptid,&result);

		if (data != NULL) free(data);
		if (weights != NULL) free(weights);
		if (bias != NULL) free(bias);
		if (signs != NULL) free(signs);
	}

	printf("Process terminating!");
	return 0;
}


/*************************** Extracts a data segment **************************/
/* Exctracts the data segment corresponding to step from matrix data. The     */
/* resulting segment is stored in windata. Chans, frames and epochs denote    */
/* the dimensions of data. The window array specifies dimensions for grid     */
/* from which to extract windata                                              */
/*                                                                            */
/* data:     doublereal array [chans,frames*epoch] (input)                    */
/* datawin:  doublereal array -                                               */
/*                   [chans,window(FRAMEWINDOW)*window(EPOCHWINDOW)] (output) */
/* chans:    int (input)                                                      */
/* frames:   int (input)                                                      */
/* epochs:   int (input)                                                      */
/* step:     int (input)                                                      */
/* window:   int array [NWINDOW] (input)                                      */

doublereal *extract(doublereal *data, doublereal *windata, int chans, int frames, int epochs, int step, int *window) {
	int i, j, nep, ofs, size;

	nep = (epochs-window[EPOCHWINDOW])/window[EPOCHSTEP] + 1;
	ofs = frames * window[EPOCHSTEP] * (step%nep);
	ofs += window[FRAMESTEP] * (step/nep);
	size = window[FRAMEWINDOW]*chans;
	
	for (i=0,j=0 ; i<window[EPOCHWINDOW] ; i++,j+=size,ofs+=frames)
		memcpy(&windata[j],&data[ofs*chans],size*sizeof(doublereal));

	return windata;
}


/*************************** Extracts data baseline ***************************/
/* Exctracts and concatenates all baseline data contained in matrix data. The */
/* resulting data block is stored in basedata. Chans, frames and epochs       */
/* denote the dimensions of data. The window array contains, amongst others,  */
/* the baseline length.                                                       */
/*                                                                            */
/* data:     doublereal array [chans,frames*epoch] (input)                    */
/* basedata: doublereal array [chans,window(BASELINE)*epoch] (output)         */
/* chans:    int (input)                                                      */
/* frames:   int (input)                                                      */
/* epochs:   int (input)                                                      */
/* window:   int array [NWINDOW] (input)                                      */

doublereal *baseline(doublereal *data, doublereal *basedata, int chans, int frames, int epochs, int *window) {
	int i, j, k, size;

	size = window[BASELINE]*chans;
	for (i=0,j=0,k=0 ; i<epochs ; i++,j+=size,k+=frames*chans)
		memcpy(&basedata[j],&data[k],size*sizeof(doublereal));

	return basedata;
}


/********************** Converts speed to CPU quantities **********************/
/* Converts speed values as represented by host info to CPU quantities in     */
/* multi-processor systems.                                                   */
/*                                                                            */
/* speed: int (input)                                                         */

int speed2proc(int speed) {
	if (speed >= 500) return speed/500;
	return 2;
}


/************************** Executes ICA in parallel **************************/
/* A serial ICA decomposition is initially performed on the baseline portions */
/* of data, all succeeding decompositions are use this result as initial      */
/* weight estimate. A number of tasks are spawned on each node according to   */
/* host info speed values. Each node is repeatedly assigned jobs until all    */
/* decompositions in the frame/epoch grid specified by window are distributed.*/
/* Weights returned by slave processes are continously sorted by projected    */
/* variance stored and to the files specified by fnames. The following        */
/* externally accessible variables are assumed initialized: extended,         */
/* extblocks, pdfsize, nsub, verbose, block, maxsteps, lrate, annealstep,     */
/* annealdeg, nochange, and momentum. If the boolean variable extended is set,*/
/* signs must be defined (i.e. not NULL)                                      */
/*                                                                            */
/* data:    double array [ncomps,frames*epoch] (input)                        */
/* weights: double array [ncomps,ncomps] (input)                              */
/* sphere:  double array [ncomps,ncomps] (input)                              */
/* eigv:    double array [ncomps,chans] (input)                               */
/* chans:   integer (input)                                                   */
/* ncomps:  integer (input)                                                   */
/* frames:  integer (input)                                                   */
/* epochs:  integer (input)                                                   */
/* window:  integer array [NWINDOW] (input)                                   */
/* bias:    double array [ncomps] (dummy) or NULL                             */
/* signs:   integer array [ncomps] (dummy) or NULL                            */
/* fnames:  char array [3] of array                                           */

void pvmica(doublereal *data, doublereal *weights, doublereal *sphere, doublereal *eigv, integer chans, integer ncomps, integer frames, integer epochs, int *window, doublereal *bias, integer *signs, char **fnames) {
	struct pvmhostinfo *hinfo;
	int                 i, j, datasize, maxep, maxfr, segs, speed, spwnd = 0, id = 0;
	int                *tids, tid, nhost, narch, nproc = 0, ntask = 0, last = 0;
	char               *name;
	FILE               *fids[3];
	integer           **srec;
	doublereal        **wrec, **brec, *basedata, *windata, *prjdata, deflr = lrate;
	tresult             result;
	tassign             assign;

/*	pvm_catchout(stdout);*/
	pvm_config(&nhost,&narch,&hinfo);
	for (i=0 ; i<nhost ; i++) nproc += speed2proc(hinfo[i].hi_speed);

	basedata = (doublereal*)malloc(ncomps*epochs*window[BASELINE]*sizeof(doublereal));
	baseline(data,basedata,ncomps,frames,epochs,window);
	runica(basedata,weights,ncomps,1,window[BASELINE]*epochs,bias,signs);
	free(basedata);

	assign.weights = weights;
	assign.chans = ncomps;
	assign.frames = window[FRAMEWINDOW];
	assign.epochs = window[EPOCHWINDOW];
	assign.bias = (int)(bias!=NULL);
	assign.signs = (int)(signs!=NULL);
	assign.extended = extended;
	assign.extblocks = extblocks;
	assign.pdfsize = pdfsize;
	assign.nsub = nsub;
	assign.verbose = verbose;
	assign.block = block;
	assign.maxsteps = maxsteps;
	assign.lrate = deflr;
	assign.annealstep = annealstep;
	assign.annealdeg = annealdeg;
	assign.nochange = nochange;
	assign.momentum = momentum;
	
	datasize = window[FRAMEWINDOW] * window[EPOCHWINDOW] * ncomps;
	windata = (doublereal*)malloc(datasize*sizeof(doublereal));
	prjdata = (doublereal*)malloc(datasize*sizeof(doublereal));
	
	maxep = (epochs-window[EPOCHWINDOW])/window[EPOCHSTEP] + 1;
	maxfr = (frames-window[FRAMEWINDOW])/window[FRAMESTEP] + 1;
	segs = maxep*maxfr;

	if (segs > nproc) {
		tids = (int*)malloc(nproc*sizeof(int));
		for (i=0 ; i<nhost ; i++) {
			name = hinfo[i].hi_name;
			speed = hinfo[i].hi_speed;
			spwnd = pvm_spawn(SPAWN_ICA,NULL,1,name,speed2proc(speed),&tids[ntask]);
			if (spwnd <= 0) error("Failed to spawn processes");
			ntask += spwnd;
		}
	}
	else {
		tids = (int*)malloc(segs*sizeof(int));
		spwnd += pvm_spawn(SPAWN_ICA,NULL,0,"",segs,tids);
		if (spwnd <= 0) error("Failed to spawn processes");
		ntask += spwnd;
	}

	for (i=0,j=0 ; i<ntask ; i++,j+=2,id++) {
		if (j >= ntask) j=1;
		assign.id   = id;
		assign.data = extract(data,windata,(int)ncomps,(int)frames,(int)epochs,(int)id,window);
		send_assign(tids[j],&assign);
	}

	wrec = (doublereal**)malloc(segs*sizeof(doublereal*));
	brec = (doublereal**)malloc(segs*sizeof(doublereal*));
	srec = (integer**)malloc(segs*sizeof(integer*));
	for (i=0 ; i<segs ; i++) {
		wrec[i] = NULL;
		brec[i] = NULL;
		srec[i] = NULL;
	}

	for (i=0 ; i<3 ; i++) fids[i] = NULL;
	if (fnames[0] != NULL) fids[0] = fopen(fnames[0],"wb");
	if (fnames[1] != NULL) fids[1] = fopen(fnames[1],"wb");
	if (fnames[2] != NULL) fids[2] = fopen(fnames[2],"wt");

	for (i=0 ; i<segs ; i++) {
		tid = receive_result(&result);
		
		printf("Received id %d\n",result.id);
		wrec[result.id] = result.weights;
		brec[result.id] = result.bias;
		srec[result.id] = result.signs;

		if (ntask < segs) {
			assign.id = id;
			assign.data = extract(data,windata,(int)ncomps,(int)frames,(int)epochs,(int)id,window);
			send_assign(tid,&assign);
			ntask++;
			id++;
		}
		else {
			send_kill(tid);
		}

		datasize = window[FRAMEWINDOW] * window[EPOCHWINDOW];
		extract(data,windata,(int)ncomps,(int)frames,(int)epochs,(int)result.id,window);
		geproj(windata,result.weights,(integer)ncomps,(integer)datasize,prjdata);
		if (eigv)
			varsort(prjdata,result.weights,sphere,&eigv[chans*(chans-ncomps)],result.bias,result.signs,(integer)ncomps,(integer)datasize,(integer)chans);
		else
			varsort(prjdata,result.weights,sphere,NULL,result.bias,result.signs,(integer)ncomps,(integer)datasize,(integer)chans);

		while (last<segs && wrec[last]!=NULL) {
			if (fids[0]!=NULL && wrec[last]!=NULL) fbc_matwrite(fids[0],chans*ncomps,wrec[last]);
			if (fids[1]!=NULL && brec[last]!=NULL) fbc_matwrite(fids[1],ncomps,brec[last]);
			if (fids[2]!=NULL && srec[last]!=NULL) iac_matwrite(fids[2],ncomps,srec[last]);
			if (wrec[last] != NULL) free(wrec[last]);
			if (brec[last] != NULL) free(brec[last]);
			if (srec[last] != NULL) free(srec[last]);
			last++;
		}
	}

	for (i=0 ; i<3 ; i++) 
		if (fids[i] != NULL) fclose(fids[i]);

	if (wrec != NULL) free(wrec);
	if (brec != NULL) free(brec);
	if (srec != NULL) free(srec);
	if (windata != NULL) free(windata);
	if (prjdata != NULL) free(prjdata);
	if (tids != NULL)    free(tids);
}
