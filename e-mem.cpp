/* ================================================================= *
 *  e-mem.cpp : Main program                                         *
 *                                                                   *
 *  E-MEM: An efficient (MUMmer-like) tool to retrieve Maximum Exact *
 *         Matches using hashing based algorithm                     *
 *                                                                   *
 *  Copyright (c) 2014, Nilesh Khiste                                *
 *  All rights reserved                                              *
 *                                                                   *
 *  This program is free software: you can redistribute it and/or    *
 *  modify it under the terms of the GNU General Public License as   *
 *  published by the Free Software Foundation, either version 3 of   *
 *  the License, or (at your option) any later version.              *
 *                                                                   *
 *  This program is distributed in the hope that it will be useful,  *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of   *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
 *  GNU General Public License for more details.                     *
 *                                                                   *
 *  You should have received a copy of the GNU General Public        *
 *  License along with this program.                                 *
 *                                                                   *
 *  This file is subject to the terms and conditions defined in the  *
 *  file 'LICENSE', which is part of this source code package.       *
 * ================================================================= */

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <unordered_map>
#include <map>
#include <vector>
#include <iterator>
#include <omp.h>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <sys/stat.h>

#include "e-mem.h"
#include "file.h"

using namespace std;
using namespace boost;

/*
 * Function builds a kmer hash for a reference sequence.
 * Input: empty refHash
 * Output: populated refHash
 */
void buildRefHash(Knode* &refHash, uint64_t totalBits, seqFileReadInfo &RefFile)
{
    uint64_t j=0;
    uint64_t currKmerPos=0, currKmer=0;
    int32_t offset=0;
    int nextKmerPosition = commonData::minMemLen - commonData::kmerSize + 2;

    vector<mapObject>::iterator it;
    it = upper_bound(RefFile.blockOfNs.begin(), RefFile.blockOfNs.end(), currKmerPos, mapObject());
    while (currKmerPos<=totalBits)
    {
        if (currKmerPos + commonData::kmerSize - 2 > totalBits)
            break;

        if(RefFile.checkKmerForNs(currKmerPos, it)){
            currKmerPos+=nextKmerPosition; // Move L-K+2 bits = 50-28+1=23 char = 46 bits
            continue;
        }

        offset = currKmerPos%DATATYPE_WIDTH;
        j=currKmerPos/DATATYPE_WIDTH; // next loc in binReads

        currKmer = RefFile.binReads[j]; // Get the unsigned 64-bit
        currKmer <<= offset; // Get base by shifting left by offset

        if (offset > DATATYPE_WIDTH-commonData::kmerSize) // Kmer split in two integers
            // Append the remaining bits from the next bin to make the full kmer
            currKmer |= ((RefFile.binReads[j+1] & global_mask_left[(commonData::kmerSize-(DATATYPE_WIDTH-offset))/2 -1])>>(DATATYPE_WIDTH-offset));
        else
            currKmer &= global_mask_left[commonData::kmerSize/2 - 1];
        /* Add kmer to the hash table */
        refHash->addKmerNode(currKmer, currKmerPos);
        currKmerPos+=nextKmerPosition; // Move L-K+2 bits = 50-28+1=23 char = 46 bits
    }
}

/*
 * Function extends the kmer match in left/right direction for
 * possible MEMs.
 * Input: currRPos : current position of matching reference Kmer
 * Input: currRPos : current position of matching query Kmer
 * Input: totalRBits : total number of bits in reference
 * Input: totalQBits : total number of bits in query
 * Input: name : reference sequence string for output
 *
 */
void helperReportMem(uint64_t &currRPos, uint64_t &currQPos, uint64_t totalRBits, uint64_t totalQBits, seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, tmpFilesInfo &arrayTmpFile, mapObject &RefNpos, mapObject &QueryNpos, vector<SeqPos> &SeqPosVec)
{
    /*
     * lRef and lQue are local variables for left extension of
     * reference and query sequence respectively. rRef and rQue
     * are their right counterparts.
     */
    uint64_t lRef=currRPos, lQue=currQPos; // Keeping lRef on currRPos-this makes offset computation simpler
    uint64_t offsetR=0,offsetQ=0;
    uint64_t rRef=currRPos+commonData::kmerSize, rQue=currQPos+commonData::kmerSize; // one character ahead of current match
    uint64_t currR=0, currQ=0;
    int i=0,j=0,mismatch=0;
    uint64_t matchSize=0;

    // Gets the real N bounds outside the k-mer
    if (!(((QueryNpos.left==0x1)?true:QueryNpos.left<=lQue) && rQue<=QueryNpos.right))
        QueryFile.getKmerLeftnRightBoundForNs(lQue, QueryNpos);

    if (!(((RefNpos.left==0x1)?true:RefNpos.left<=lRef) && rRef<=RefNpos.right))
        RefFile.getKmerLeftnRightBoundForNs(lRef, RefNpos);

    if (RefNpos.right-((RefNpos.left==0x1)?0:RefNpos.left)+2 < static_cast<uint64_t>(commonData::minMemLen))
        return;

    if (QueryNpos.right-((QueryNpos.left==0x1)?0:QueryNpos.left)+2 < static_cast<uint64_t>(commonData::minMemLen))
        return;

    //match towards left
    while (lRef && lQue && ((QueryNpos.left==0x1)?true:QueryNpos.left<=lQue) && ((RefNpos.left==0x1)?true:RefNpos.left<=lRef))
    {
        if (!mismatch)
        {
            offsetR=(lRef)%DATATYPE_WIDTH;
            i=(lRef)/DATATYPE_WIDTH;
            offsetQ=(lQue)%DATATYPE_WIDTH;
            j=(lQue)/DATATYPE_WIDTH;
            // Get the matchsize from reference or query that is left in the 64 bit
            if (offsetR > offsetQ)
                matchSize = offsetQ;
            else
                matchSize = offsetR;

            if (!matchSize)
                matchSize=2;

            if ((QueryNpos.left!=0x1) && (matchSize > lQue-QueryNpos.left))
                matchSize = lQue-QueryNpos.left;

            if ((RefNpos.left!=0x1) && (matchSize > lRef-RefNpos.left))
                matchSize = lRef-RefNpos.left;

            if (!matchSize)
                break;

            /*
             * There will never be case with offset=0 and i=0 because
             * i=0 happens only when lRef=0 and in that case we do not
             * enter this loop.
             */
            // Get the 64-bit representation
            currR = RefFile.binReads[offsetR?i:i-1];
            if (offsetR)
                currR >>= DATATYPE_WIDTH-offsetR;
            currQ = QueryFile.binReads[offsetQ?j:j-1];
            if (offsetQ)
                currQ >>= DATATYPE_WIDTH-offsetQ;
        }

        // if no exact match between the rest of query and ref 64 bit offset then allow for mismatch
        if((currR & global_mask_right[matchSize/2 - 1]) != (currQ &  global_mask_right[matchSize/2 - 1])) {
            if (matchSize==2)
                break;

            mismatch=1;
            matchSize/=2;
            if (matchSize%2)
                matchSize+=1;
        // if there is an exact match then extend previous lRef and lQue to the left
        }else {
            lRef-=matchSize;
            lQue-=matchSize;
            if (mismatch) { // if already had a mismatch
                if (matchSize==2) // and there is no more to extend, break
                    break;
                currR >>= matchSize;
                currQ >>= matchSize;
            }
        }
    } // loop until extension reached to the left

    if (totalRBits-lRef+2 < static_cast<uint64_t>(commonData::minMemLen))
        return;

    if (totalQBits-lQue+2 < static_cast<uint64_t>(commonData::minMemLen))
        return;

    //match towards right - same principal as the left
    mismatch=0;
    while ((rRef <= totalRBits) && (rQue <= totalQBits) && (rRef <= RefNpos.right) && (rQue <= QueryNpos.right))
    {
        if (!mismatch)
        {
            offsetR=rRef%DATATYPE_WIDTH;
            i=rRef/DATATYPE_WIDTH;
            offsetQ=rQue%DATATYPE_WIDTH;
            j=rQue/DATATYPE_WIDTH;

            if (offsetR > offsetQ)
                matchSize = DATATYPE_WIDTH-offsetR;
            else
                matchSize = DATATYPE_WIDTH-offsetQ;

            if (rRef+matchSize > totalRBits)
                matchSize = totalRBits-rRef;

            if (rQue+matchSize > totalQBits)
                matchSize = totalQBits-rQue;

            if (rQue+matchSize > QueryNpos.right)
                matchSize = QueryNpos.right-rQue;

            if (rRef+matchSize > RefNpos.right)
                matchSize = RefNpos.right-rRef;

            if(!matchSize)
                matchSize=2;



            currR = RefFile.binReads[i];
            currR <<= offsetR;
            currQ = QueryFile.binReads[j];
            currQ <<= offsetQ;
        }

        if((currR & global_mask_left[matchSize/2 - 1]) != (currQ &  global_mask_left[matchSize/2 - 1])) {
            if (matchSize==2){
                rRef-=2;
                rQue-=2;
                break;
            }
            /* if current rRef/rQue plus matchSize smaller than minMEMLength, then simply return.
             * Note that one less character is compared due to a mismatch
             */
            if (rRef+matchSize-lRef < static_cast<uint64_t>(commonData::minMemLen))
            {
                return;
            }

            mismatch=1;
            matchSize/=2;
            if (matchSize%2)
                matchSize+=1;
        }else {
            if (mismatch) {
                if (matchSize==2)
                    break;
            }
            if ((rRef == totalRBits) || (rQue == totalQBits))
                break;

            currR <<= matchSize;
            currQ <<= matchSize;
            rRef+=matchSize;
            rQue+=matchSize;
        }
    } // loop until extension reached to the right

    /* Adjust rRef and rQue locations */

    if (rRef > RefNpos.right){
        rQue-=(rRef-RefNpos.right);
        rRef=RefNpos.right;
    }
    if (rQue > QueryNpos.right){
        rRef-=(rQue-QueryNpos.right);
        rQue=QueryNpos.right;
    }

    if (rRef > totalRBits){
        rQue-=(rRef-totalRBits);
        rRef=totalRBits;
    }
    if (rQue > totalQBits){
        rRef-=(rQue-totalQBits);
        rQue=totalQBits;
    }

    /* Ignore reverse complements of the same ORF, i.e. where matching prefix/suffix of reference/query
     * Also excludes N mismatches
     */
    if ((lRef?lRef!=RefNpos.left:!RefNpos.left) && (rQue?rQue!=QueryNpos.right:!QueryNpos.right)){
        if ((rRef?rRef!=RefNpos.right:!RefNpos.right) && (lQue?lQue!=QueryNpos.left:!QueryNpos.left)) {
            for (vector<SeqPos>::iterator it=SeqPosVec.begin();it!=SeqPosVec.end();++it) {
                if ((*it).lR == lRef && (*it).rR == rRef && (*it).lQ == lQue && (*it).rQ == rQue) {
                    return;
                }
            }
            arrayTmpFile.writeMemInTmpFiles(lRef, rRef, lQue, rQue, RefNpos, QueryNpos, QueryFile, RefFile);
            SeqPos seqPos(lRef, rRef, lQue, rQue);
            SeqPosVec.push_back(seqPos);
        }
    }
}

void reportMEM(Knode* &refHash, uint64_t totalBases, uint64_t totalQBases, seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, tmpFilesInfo &arrayTmpFile)
{
    uint64_t totalQBits = CHARS2BITS(totalQBases); // convert char position to bit position of QueryFile.totalBases-1
    uint32_t copyBits=0;
    #pragma omp parallel num_threads(commonData::numThreads)
    {
        uint64_t currKmer=0, j=0;
        int32_t offset=0;
        uint32_t first=1;
        int kmerWithNs=0;
        mapObject QueryNpos, RefNpos;
        vector<SeqPos> SeqPosVec;
        vector<mapObject>::iterator it;
        it = upper_bound(QueryFile.blockOfNs.begin(), QueryFile.blockOfNs.end(), 0, mapObject());

        #pragma omp single
        {
        /*
        * Number of copy bits during query kmer processing depends on kmer size.
        */
        if (DATATYPE_WIDTH-commonData::kmerSize > 32 )
            copyBits=32; //16 characters
        else if (DATATYPE_WIDTH-commonData::kmerSize > 16) //kmerSize = 40 so copyBits = 16
            copyBits=16; //8 characters
        else
            copyBits=8; //4 characters

        /* If copyBits more than 8, the for loop parallelisation will give
        * incorrect results - miss some Mems
        */
        if(commonData::numThreads > 1)
            copyBits=8; //4 characters // if more than one thread, then copyBits = 8
        }

        // parallelise search of query kmer in reference
        #pragma omp for
        for (uint64_t currKmerPos=0; currKmerPos<=totalQBits; currKmerPos+=2)
        {
            if ((currKmerPos + commonData::kmerSize - 2) > totalQBits)
                continue;

            if(QueryFile.checkKmerForNs(currKmerPos, it)){
                kmerWithNs=1;
            }

            j=currKmerPos/DATATYPE_WIDTH;// current location in binReads
            offset = currKmerPos%DATATYPE_WIDTH;
            if(first || !offset){
                currKmer = QueryFile.binReads[j];
                currKmer <<= offset;
                if(offset > DATATYPE_WIDTH-commonData::kmerSize)
                  currKmer |= ((QueryFile.binReads[j+1] & global_mask_left[offset/2-1])>>(DATATYPE_WIDTH-offset));
                first=0;
            }else
                currKmer <<= 2;

            if(offset  && !(offset % copyBits))
              currKmer |= ((QueryFile.binReads[j+1] & global_mask_left[offset/2-1])>>(DATATYPE_WIDTH-offset));

            if (kmerWithNs){
                /* Do not process this Kmer, Ns in it */
                kmerWithNs=0;
                continue;
            }
            /* Find the K-mer in the refHash */
            uint64_t *dataPtr=NULL;
            if (refHash->findKmer(currKmer & global_mask_left[commonData::kmerSize/2 - 1], dataPtr)) // dataPtr is position of the kmer in reference
            {
                // We have a match
                for (uint64_t n=1; n<=dataPtr[0]; n++) { // currKmerPos is position of kmer in query
                    helperReportMem(dataPtr[n], currKmerPos, CHARS2BITS(totalBases), CHARS2BITS(totalQBases), RefFile,QueryFile, arrayTmpFile, RefNpos, QueryNpos, SeqPosVec);
                }
            }
        }
    }
}

void searchQuery(Knode* &refHash, seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, tmpFilesInfo &arrayTmpFile)
{
    QueryFile.resetCurrPos();
    reportMEM(refHash, RefFile.totalBases-1, QueryFile.totalBases-1, RefFile, QueryFile, arrayTmpFile);
    // QueryFile.setCurrPos();
    // QueryFile.clearMapForNs();
    // QueryFile.clearTmpString();
}

void processReference(seqFileReadInfo &RefFile, seqFileReadInfo &QueryFile, tmpFilesInfo &arrayTmpFile)
{
    uint64_t numberOfKmers=0,n=0;
    int hashTableSizeIndex=0;
    Knode *refHash;

    numberOfKmers = ceil((RefFile.totalBases-commonData::kmerSize/2+1)/((commonData::minMemLen/2-commonData::kmerSize/2 + 1)) + 1);
    //cout << "number of kmers: " << numberOfKmers << endl;

    /* Set the size of the hash table to the numberofKmers. */
    for (n=0; n<450; ++n)
    {
        if (hashTableSize[n] > 1.75*numberOfKmers)
        {
            hashTableSizeIndex = n;
            break;
        }
    }
    // cout << "Hashtable size: " << hashTableSize[hashTableSizeIndex] << endl;

    Knode::currHashTabSize = hashTableSize[hashTableSizeIndex];  //Store the size of the hash table.
    if (hashTableSizeIndex)
        Knode::prevHashTabSize = hashTableSize[hashTableSizeIndex-1];
    else
        Knode::prevHashTabSize = 3;

    /* Create the refHash for K-mers. */
    cout << "allocate hashtable" << endl;
    refHash = new Knode[Knode::currHashTabSize];
    cout << "build hashtable" << endl;
    buildRefHash(refHash, CHARS2BITS(RefFile.totalBases-1), RefFile);
    cout << "Query search" << endl;
    searchQuery(refHash, RefFile, QueryFile, arrayTmpFile);

    delete [] refHash;
}

bool is_numeric(const string &str)
{
    return all_of(str.begin(), str.end(), ::isdigit);
}

void checkCommandLineOptions(uint32_t &options)
{
    if (!IS_FASTAU_DEF(options)){
        if (!IS_FASTA1_DEF(options) || !IS_FASTA2_DEF(options)){
            cout << "ERROR: both f1 forward and f2 reverse fasta files or fu unpaired fasta file must be passed!" << endl;
            exit(EXIT_FAILURE);
        }
    }

    if (!IS_OUT_FILE_DEF(options)) {
        cout << "ERROR: output file must be passed!" << endl;
        exit(EXIT_FAILURE);
    }

    if (IS_SPLIT_SIZE_DEF(options)){
        if (commonData::d <= 0){
            cout << "ERROR: -d cannot be less than or equal to zero!" << endl;
            exit(EXIT_FAILURE);
        }
    }

    if (IS_NUM_THREADS_DEF(options)){
        if (commonData::numThreads <= 0){
            cout << "ERROR: -t cannot be less than or equal to zero!" << endl;
            exit(EXIT_FAILURE);
        }
    }

    if (IS_LENGTH_DEF(options)){
        if (commonData::minMemLen <= 2){
            cout << "ERROR: -l cannot be less than or equal to one!" << endl;
            exit(EXIT_FAILURE);
        }
        if ((commonData::minMemLen - commonData::kmerSize) < 0)
        {
            if (IS_KMER_SIZE_DEF(options)){
                cout << "ERROR: kmer-size cannot be larger than min mem length!" << endl;
                exit( EXIT_FAILURE );
            }else {
                commonData::kmerSize = commonData::minMemLen;
            }
        }
    }

    if (IS_KMER_SIZE_DEF(options)){
        if (commonData::kmerSize <= 0){
            cout << "ERROR: -k cannot be less than or equal to zero!" << endl;
            exit(EXIT_FAILURE);
        }

        if (commonData::kmerSize > CHARS2BITS(28))
        {
            cout << "ERROR: kmer-size cannot be greater than 28" << endl;
            exit( EXIT_FAILURE );
        }
    }
}

void print_help_msg()
{
    cout <<  endl;
    cout << "PAL-MEM Version 1.0.0, Jan. 31, 2020" << endl;
    cout << "Adapted from E-MEM Version 1.0.2, Dec. 12, 2017, by Nilesh Khiste and Lucian Ilie" << endl;
    cout <<  endl;
    cout << "E-MEM finds and outputs the position and length of all maximal" << endl;
    cout << "exact matches (MEMs) between <query-file> and <reference-file>" << endl;
    cout << endl;
    cout << "Usage: e-mem [options]" << endl;
    cout << endl;
    cout << "Options:" << endl;
    cout << "-f1\t<filename>\t" << "fasta file with forward paired-end reads" << endl;
    cout << "-f2\t<filename>\t" << "fasta file with reverse paired-end reads" << endl;
    cout << "-fu\t<filename>\t" << "fasta file with unpaired reads" << endl;
    cout << "-o\t<filename>\t" << "output file" << endl;
    cout << "-l\t" << "set the minimum length of a match. The default length" << endl;
    cout << "  \tis 50" << endl;
    cout << "-c\t" << "report the query-position of a reverse complement match" << endl;
    cout << "  \trelative to the original query sequence" << endl;
    cout << "-F\t" << "force 4 column output format regardless of the number of" << endl;
    cout << "  \treference sequence input" << endl;
    cout << "-L\t" << "show the length of the query sequences on the header line" << endl;
    cout << "-d\t" << "set the split size. The default value is 1" << endl;
    cout << "-t\t" << "number of threads. The default is 1 thread" << endl;
    cout << "-h\t" << "show possible options" << endl;
}

int main (int argc, char *argv[])
{
    int32_t i=0, n=1;
    uint32_t options=0;
    seqFileReadInfo QueryFile, RefFile;
    outFileReadInfo OutFile;
    string fasta1, fasta2, fastaU, outFilename;

    // Check Arguments
    if (argc==1 || argc==2){
        print_help_msg();
        exit(EXIT_SUCCESS);
    }

    while(argv[n]) {

        if(boost::equals(argv[n], "-f1")){
            if (IS_FASTA1_DEF(options)){
                cout << "ERROR: f1 argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            if (IS_FASTAU_DEF(options)){
                cout << "ERROR: the f1 and f2 arguments must be used together or fu argument used alone!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_FASTA1(options);
            fasta1 = argv[n+1];
            n+=2;
        }else if(boost::equals(argv[n], "-f2")){
            if (IS_FASTA2_DEF(options)){
                cout << "ERROR: f2 argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            if (IS_FASTAU_DEF(options)){
                cout << "ERROR: the f1 and f2 arguments must be used together or fu argument used alone!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_FASTA2(options);
            fasta2 = argv[n+1];
            n+=2;
        }else if(boost::equals(argv[n], "-fu")) {
            if (IS_FASTAU_DEF(options)) {
                cout << "ERROR: fu argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            if (IS_FASTA1_DEF(options) || IS_FASTA2_DEF(options)) {
                cout << "ERROR: the f1 and f2 arguments must be used together or fu argument used alone!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_FASTAU(options);
            fastaU = argv[n + 1];
            n += 2;
        } else if(boost::equals(argv[n], "-o")){
            if (IS_OUT_FILE_DEF(options)) {
                cout << "ERROR: output file argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_OUT_FILE(options);
            outFilename = argv[n+1];
            n+=2;
        }else if(boost::equals(argv[n],"-l")){
            if (IS_LENGTH_DEF(options)) {
                cout << "ERROR: Length argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_LENGTH(options);
            if (!argv[n+1] || !is_numeric(argv[n+1])){
                cout << "ERROR: Invalid value for -l option!" << endl;
                exit(EXIT_FAILURE);
            }
            commonData::minMemLen = 2*std::stoi(argv[n+1]);
            n+=2;
        }else if (boost::equals(argv[n],"-d")){
            if (IS_SPLIT_SIZE_DEF(options)) {
                cout << "ERROR: Split size argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            if (!argv[n+1] || !is_numeric(argv[n+1])){
                cout << "ERROR: Invalid value for -d option!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_SPLIT_SIZE(options);
            commonData::d = std::stoi(argv[n+1]);
            n+=2;
        }else if (boost::equals(argv[n],"-t")){
            if (IS_NUM_THREADS_DEF(options)) {
                cout << "ERROR: Number of threads argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            if (!argv[n+1] || !is_numeric(argv[n+1])){
                cout << "ERROR: Invalid value for -t option!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_NUM_THREADS(options);
            commonData::numThreads = std::stoi(argv[n+1]);
            n+=2;
        }else if (boost::equals(argv[n],"-k")){
            if (IS_KMER_SIZE_DEF(options)) {
                cout << "ERROR: Kmer size argument passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            if (!argv[n+1] || !is_numeric(argv[n+1])){
                cout << "ERROR: Invalid value for -k option!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_KMER_SIZE(options);
            commonData::kmerSize = 2*std::stoi(argv[n+1]);
            n+=2;
        }else if (argv[n][0] != '-'){
            cout << "ERROR: option must start with '-'!" << endl;
            exit(EXIT_FAILURE);
        }else if (boost::equals(argv[n],"-c")){
            if (IS_RELREV_QUEPOS_DEF(options)) {
                cout << "ERROR: option -c passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_RELREV_QUEPOS(options);
            n+=1;
            commonData::relQueryPos = 1;
        }else if (boost::equals(argv[n],"-F")){
            if (IS_FCOL_OUTPUT_DEF(options)) {
                cout << "ERROR: option -F passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_FCOL_OUTPUT(options);
            n+=1;
            commonData::fourColOutput = 1;
        }else if (boost::equals(argv[n],"-L")){
            if (IS_LEN_IN_HEADER_DEF(options)) {
                cout << "ERROR: option -L passed multiple times!" << endl;
                exit(EXIT_FAILURE);
            }
            SET_LEN_IN_HEADER(options);
            n+=1;
            commonData::lenInHeader = 1;
        }else if (boost::equals(argv[n],"-mum")){
            /* Ignore this option and continue */
            n+=1;
        }else if (boost::equals(argv[n],"-mumcand")){
            /* Ignore this option and continue */
            n+=1;
        }else if (boost::equals(argv[n],"-mumreference")){
            /* Ignore this option and continue */
            n+=1;
        }else if (boost::equals(argv[n],"-maxmatch")){
            /* Ignore this option and continue */
            n+=1;
        }else if (boost::equals(argv[n],"-s")){
            /* Ignore this option and continue */
            n+=1;
        }else if (boost::equals(argv[n],"-h")){
            print_help_msg();
            exit(EXIT_SUCCESS);
        }else {
            cout << "ERROR: Invalid option." << endl << flush;
            print_help_msg();
            exit( EXIT_FAILURE );
        }
    }

    checkCommandLineOptions(options);

    // Open tmp files
    cout << "Open tmp files" << endl;
    sprintf(commonData::nucmer_path, "%s/%d_tmp", getenv("NUCMER_E_MEM_OUTPUT_DIRPATH")?getenv("NUCMER_E_MEM_OUTPUT_DIRPATH"):".",getpid());

    tmpFilesInfo arrayTmpFile(NUM_TMP_FILES+2);
    arrayTmpFile.openFiles(ios::out|ios::binary, NUM_TMP_FILES+2);

    cout << "Open fasta files" << endl;
    if (IS_FASTA1_DEF(options) && IS_FASTA2_DEF(options)){
        vector<string> filenames = {fasta1, fasta2};
        QueryFile.setFiles(filenames);
        RefFile.setFiles(filenames);
    } else if (IS_FASTAU_DEF(options)) {
        vector<string> filenames = {fastaU};
        QueryFile.setFiles(filenames);
        RefFile.setFiles(filenames);
    }

    cout << "Generate reverse complement" << endl;
    QueryFile.generateRevComplement();
    //RefFile.generateRevComplement(0); // This routine also computes size and num sequences
    RefFile.setSize(QueryFile.getSize());
    RefFile.setNumSequences(QueryFile.getNumSequences());

    /* Only reverse complement matches */
    QueryFile.setReverseFile();

    arrayTmpFile.setNumMemsInFile(QueryFile.allocBinArray(0), QueryFile.getNumSequences());
    //cout << "Size: " << QueryFile.getSize() << endl;
    RefFile.allocBinArray(1);
    //cout << "Size: " << RefFile.getSize() << endl;
    RefFile.clearFileFlag();

    QueryFile.clearFileFlag();
    cout << "Encode Query sequences" << endl;
    if (QueryFile.readChunks())
    {
        for (i=0; i<commonData::d; i++) {
            cout << "Encode Ref sequences" << endl;
            if(RefFile.readChunks()){ // Encode sequence as 2-bits in RefFile object
                cout << "Build Ref hashtable" << endl;
                processReference(RefFile, QueryFile, arrayTmpFile); // Build hashtable, query hashtable, find ls, and write temp files
                RefFile.setCurrPos();
                RefFile.clearMapForNs(); // clear block of Ns from memory
            }
            else
                break;
        }

        /*
         * Process MemExt list and write to file
         */
        cout << "If MemExtVec, merge to make temp files" << endl;
        arrayTmpFile.mergeMemExtVector();
    }

    /*
     * Free up the allocated arrays
     */
    arrayTmpFile.closeFiles(NUM_TMP_FILES);
    QueryFile.closeFile();
    QueryFile.destroy();

    /*
     * Populate sequence information in vectors. Use this to get MEM
     * Populate sequence information in vectors. Use this to get MEM
     * positions relative to the original sequences.
     */
    cout << "Create seqData" << endl;
    vector<seqData> refSeqInfo;
    // vector<seqData> querySeqInfo;
    refSeqInfo.reserve(RefFile.getNumSequences());
    // querySeqInfo.reserve(QueryFile.getNumSequences());
    RefFile.generateSeqPos(refSeqInfo);
    //QueryFile.generateSeqPos(querySeqInfo);

    cout << "Get inverted repeats" << endl;
    vector<posData> invertedRepeatInfo;
    arrayTmpFile.getInvertedRepeats(RefFile, refSeqInfo, invertedRepeatInfo);

    cout << "Write inverted repeats" << endl;
    arrayTmpFile.writeInvertedRepeats(RefFile, refSeqInfo, invertedRepeatInfo);

    cout << "Write to outfile" << endl;
    OutFile.setFile(outFilename);
    if (IS_FASTA1_DEF(options) && IS_FASTA2_DEF(options)){
        vector<string> filenames = {fasta1, fasta2};
        OutFile.removeSeq(refSeqInfo, filenames);
    } else if (IS_FASTAU_DEF(options)) {
        vector<string> filenames = {fastaU};
        OutFile.removeSeq(refSeqInfo, filenames);
    }

    //arrayTmpFile.removeTmp();

    RefFile.closeFile();
    OutFile.closeFile();

    RefFile.destroy();

    fflush(0);
    return 0;
}