/* ================================================================= *
 *  file.h : Header file with supporting class definitions           *
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
using namespace std;

#define CHARS2BITS(x) 		(2*(x)) //convert char position to bit position
#define DATATYPE_WIDTH          64 	// number of bits
#define RANDOM_SEQ_SIZE         10
#define NUM_TMP_FILES           24
#define HAMMING_DISTANCE_LIM    1

class commonData {
  public:
    static int32_t minMemLen;
    static int32_t d;
    static int32_t numThreads;
    static int32_t kmerSize;
    static int32_t fourColOutput;
    static int32_t lenInHeader;
    static int32_t relQueryPos;
    static char nucmer_path[256];
};

int32_t commonData::minMemLen=100; // 2 bit representation=50
int32_t commonData::d=1;
int32_t commonData::numThreads=1;
int32_t commonData::kmerSize=56; //2 bit representation = 28
int32_t commonData::fourColOutput=0;
int32_t commonData::lenInHeader=0;
int32_t commonData::relQueryPos=0;
char commonData::nucmer_path[256]={'\0'};


class seqData {
  public:
    uint64_t start;
    uint64_t end;
    std::string seq;
    int keep=1;
    seqData()
    {
        start=0;
        end=0;
    };

    bool operator ()(const seqData &obj1, const seqData &obj2)
    {
      return (obj2.start>obj1.end?true:false);
    }
};

class posData {
public:
    uint64_t L1Bound;
    uint64_t R1Bound;
    uint64_t L2Bound;
    uint64_t R2Bound;
    uint64_t palLength;
    std::string seq;

    posData() {}

    bool operator () (const posData &obj1, const posData &obj2)
    {
        if (obj1.L1Bound < obj2.L1Bound) {
            return true;
        } else if (obj1.L1Bound > obj2.L1Bound) {
            return false;
        } else {
            if (obj1.L2Bound < obj2.L2Bound)
                return true;
            else
                return false;
        }
    }
};

class mapObject {
    public:
        uint64_t left;
        uint64_t right;
        mapObject() {
          left=0;
          right=0;
        }

        mapObject(uint64_t x, uint64_t y) {
          left=x;
          right=y;
        }

        bool operator()(const uint64_t &x, const mapObject &y)
        {
          return x < y.left;
        }
};


class seqFileReadInfo {
    fstream file1, file2;
    int numSeqFiles;
    uint64_t size;
    string strTmp, strName;
    uint64_t binReadSize;
    uint64_t binReadsLocation;
    uint64_t currPos;
    uint64_t numSequences;

        string& randomStr()
        {
             static string str("NNNNNNNNNN");
             return str;
        }

        void processTmpString(uint64_t &sz, uint64_t &blockNCount)
        {
            string line = strTmp;
            strTmp.clear();
            totalBases=0;
            binReadsLocation=0;
            processInput(line, sz, blockNCount);
        }

      /*
       * Function converts a character sequence into an array of integers.
       * Input: character string
       * Output: array of integers, total number of bases
       */
    void processInput(string &str, uint64_t &sz, uint64_t &blockNCount)
    {
        int chooseLetter=0;

        /* Processing the sequences by encoding the base pairs into 2 bits. */
        for ( std::string::iterator it=str.begin(); it!=str.end(); ++it){
            if (totalBases == sz){
              strTmp += *it;
              continue;
            }else if (totalBases >= size) {
              strTmp += *it;
            }
            switch(*it)
                {
                    case 'A':
                    case 'a':
                    binReads[binReadsLocation] <<= 2; // shift left by 2 bits
                    if (blockNCount){
                        blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                        blockNCount=0;
                    }
                    break;
                    case 'C':
                    case 'c':
                        binReads[binReadsLocation] <<= 2;
                        binReads[binReadsLocation] |= 1;
                    if (blockNCount){
                        blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                        blockNCount=0;
                    }
                    break;
                    case 'G':
                    case 'g':
                        binReads[binReadsLocation] <<= 2;
                        binReads[binReadsLocation] |= 2;
                    if (blockNCount){
                        blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                        blockNCount=0;
                    }
                    break;
                    case 'T':
                    case 't':
                        binReads[binReadsLocation] <<= 2;
                        binReads[binReadsLocation] |= 3;
                    if (blockNCount){
                        blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                        blockNCount=0;
                    }
                    break;
                    default:
                    if(!blockNCount)
                        blockNCount=totalBases+1;
                        chooseLetter = rand() % 4;
                    if (chooseLetter == 0)
                        binReads[binReadsLocation] <<= 2;
                    else if (chooseLetter == 1)
                    {
                        binReads[binReadsLocation] <<= 2;
                        binReads[binReadsLocation] |= 1;
                    }
                    else if (chooseLetter == 2)
                    {
                        binReads[binReadsLocation] <<= 2;
                        binReads[binReadsLocation] |= 2;
                    }
                    else
                    {
                        binReads[binReadsLocation] <<= 2;
                        binReads[binReadsLocation] |= 3;
                    }
                }
            totalBases++;
            if ((totalBases%32)==0){
                binReadsLocation++;
            }
        }
    }

    public:
      uint64_t *binReads;
      uint64_t totalBases;
      std::vector <mapObject> blockOfNs;

      seqFileReadInfo() {
          size=0;
          currPos=0;
          binReadSize=0;
          binReadsLocation=0;
          numSequences=0;
          totalBases=0;
          numSeqFiles=0;
      }

      void setFiles(vector<string> &filenames)
      {
          numSeqFiles=filenames.size();
          if (numSeqFiles >= 1){
              file1.open(filenames[0], ios::in);
              if(!file1.is_open()) {
                  cout << "ERROR: unable to open " << filenames[0] << " file" << endl;
                  exit( EXIT_FAILURE );
              }
              if (numSeqFiles == 2){
                  file2.open(filenames[1], ios::in);
                  if(!file2.is_open()) {
                      cout << "ERROR: unable to open " << filenames[1] << " file" << endl;
                      exit( EXIT_FAILURE );
                  }
              }
          }
      }

      void setNumSequences(uint64_t numSeq) {
          numSequences = numSeq;
      }

      void setSize(uint64_t sizeSet) {
          size = sizeSet;
      }

      uint64_t &getNumSequences() {
          return numSequences;
      }

      uint64_t &getSize() {
          return size;
      }

      void openFile(string s, fstream &file){
          file.open(s, ios::in);
          if(!file.is_open()) {
              cout << "ERROR: unable to open "<< s << " file" << endl;
              exit( EXIT_FAILURE );
          }
      }

      void setReverseFile() {
          char buffer[256];
          memset(buffer,0,256);
          sprintf(buffer, "%s/revComp", commonData::nucmer_path);

          if (numSeqFiles >= 1){
              file1.close();
              if (numSeqFiles == 2){
                  file2.close();
              }
          }
          openFile(buffer, file1);
          numSeqFiles = 1;
      }

      void closeFile() {
          if (numSeqFiles >= 1){
              file1.close();
              if (numSeqFiles == 2){
                  file2.close();
              }
          }
      }

      void destroy() {
          currPos=0;
          binReadSize=0;
          binReadsLocation=0;
          totalBases=0;
          strName.clear();
          strTmp.clear();
          clearMapForNs();
          delete [] binReads;
      }

      void clearFileFlag()
      {
          if (numSeqFiles >= 1){
              file1.clear();
              file1.seekg(0, ios::beg);
              if (numSeqFiles == 2){
                  file2.clear();
                  file2.seekg(0, ios::beg);
              }
          }
      }

      uint64_t allocBinArray(int split)
      {
          if (split)
              size = size/commonData::d;
          binReadSize = floor((size+numSequences*RANDOM_SEQ_SIZE+commonData::d)/32+4);
          binReads = new uint64_t[binReadSize+1];
          memset(binReads, 0, sizeof(uint64_t)*(binReadSize+1));
          if (numSequences)
              blockOfNs.reserve(numSequences+10);
          return size;
      }

      void clearMapForNs()
      {
          blockOfNs.clear();
      }

      void clearTmpString()
      {
          strTmp.clear();
          strName.clear();
          clearMapForNs();
      }

      void getKmerLeftnRightBoundForNs(uint64_t &currKmerPos, mapObject &bounds)
      {
          uint64_t right=0;
          /*
           * Since we do all computation with bits, all our
           * positions are even. Here I return 1 (odd position),
           * an indication of no Ns towards left
           */

          if (!blockOfNs.size()){
              bounds.left=0x1;
              bounds.right=CHARS2BITS(totalBases-1);
              return;
          }

          vector<mapObject>::iterator it;
          it=upper_bound(blockOfNs.begin(), blockOfNs.end(), currKmerPos, mapObject());
          /* No N block beyond this point */
          if (it == blockOfNs.end())
              right = CHARS2BITS(totalBases-1);
          else
              right = (*it).left-2;

          /* This function never gets a position which is N */
          if (!currKmerPos || it==blockOfNs.begin()){
              bounds.left=0x1;
              bounds.right=right;
              return;
          }

          --it;

          bounds.left=(*it).right+2;
          bounds.right=right;
          return;
      }

      bool checkKmerForNs(uint64_t &currKmerPos, vector<mapObject>::iterator &it)
      {
          if (!blockOfNs.size())
              return false;

          while(it != blockOfNs.end())
          {
              if ((*it).left>currKmerPos)
                  break;
              else
                  ++it;
          }


          /* No N block beyond this point */
          if (it == blockOfNs.end()){
              --it;
              /* Current position within N block */
              if (((*it).left <=currKmerPos) && (currKmerPos <= (*it).right)){
                  ++it;
                  return true;
              }else{
                  ++it;
                  return false;
              }
          }

          if ((*it).left > (currKmerPos+commonData::kmerSize-2)){
              if (it != blockOfNs.begin()){
                  --it;
                  if ((*it).right < currKmerPos){
                      ++it;
                      return false;
                  }else {
                      ++it;
                      return true;
                  }
              }else
                  return false;
          }else {
              return true;
          }
      }


      void setCurrPos() {
          currPos+=size;;
      }

      uint64_t getCurrPos() {
          return currPos;
      }

      void resetCurrPos() {
          currPos=0;
      }

      bool readChunks()
      {
          string line;
          uint64_t blockNCount=0;
          int minSize = commonData::minMemLen/2-1;
          uint64_t sz=size+minSize;
          /* Process anything remaining from the last iteration */
          processTmpString(sz, blockNCount);

          if (numSeqFiles >= 1){
              while(getline( file1, line ).good()){
                  if(line[0] == '>' || (totalBases == sz)){
                      if( !strName.empty()){
                         if(line[0] != '>') {
                              processInput(line, sz, blockNCount);
                          }else {
                              processInput(randomStr(), sz, blockNCount);
                          }
                          if (totalBases == sz) {
                              if ((totalBases%32)!=0)
                              {
                                  uint64_t offset = CHARS2BITS(totalBases)%DATATYPE_WIDTH;
                                  binReads[binReadsLocation] <<= (DATATYPE_WIDTH-offset);
                                  binReadsLocation++;
                                  binReads[binReadsLocation]=0;
                              }
                              if (blockNCount){
                                  blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                                  blockNCount=0;
                              }
                              return true;
                          }
                      }
                      if( !line.empty() ){
                          strName = line.substr(1);
                      }
                  } else if( !strName.empty() ){
                      processInput(line, sz, blockNCount);
                  }
              }
              if (numSeqFiles == 2) {
                  while(getline( file2, line ).good()){
                      if(line[0] == '>' || (totalBases == sz)){
                          if( !strName.empty()){
                             if(line[0] != '>') {
                                  processInput(line, sz, blockNCount);
                              }else {
                                  processInput(randomStr(), sz, blockNCount);
                              }
                              if (totalBases == sz) {
                                  if ((totalBases%32)!=0)
                                  {
                                      uint64_t offset = CHARS2BITS(totalBases)%DATATYPE_WIDTH;
                                      binReads[binReadsLocation] <<= (DATATYPE_WIDTH-offset);
                                      binReadsLocation++;
                                      binReads[binReadsLocation]=0;
                                  }
                                  if (blockNCount){
                                      blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                                      blockNCount=0;
                                  }
                                  return true;
                              }
                          }
                          if( !line.empty() ){
                              strName = line.substr(1);
                          }
                      } else if( !strName.empty() ){
                          processInput(line, sz, blockNCount);
                      }
                  }
              }
          }
          if( !strName.empty() ){
              if ((totalBases%32)!=0)
              {
                  uint64_t offset = CHARS2BITS(totalBases)%DATATYPE_WIDTH;
                  binReads[binReadsLocation] <<= (DATATYPE_WIDTH-offset);
                  binReadsLocation++;
                  binReads[binReadsLocation]=0;
              }
              if (blockNCount){
                  blockOfNs.push_back(mapObject(CHARS2BITS(blockNCount-1), CHARS2BITS(totalBases-1)));
                  blockNCount=0;
              }
              if (!strTmp.size())
                  strName.clear();
              return true;
          }
          return false;
      }

      void flipCharacter(char &in, char &out)
      {
          switch(in)
          {
              case 'A':
              case 'a':
                  out='T';
                  break;
              case 'C':
              case 'c':
                  out='G';
                  break;
              case 'G':
              case 'g':
                  out='C';
                  break;
              case 'T':
              case 't':
                  out='A';
                  break;
              default:
                  out=in;
          }
      }

      void flipNswap(string &content)
      {
          string::iterator itBeg = content.begin();
          string::iterator itEnd = --content.end();
          char beg=0, end=0;
          uint64_t d=0;
          while ((d=distance(itBeg,itEnd)))
          {
              flipCharacter(*itBeg, end);
              flipCharacter(*itEnd, beg);
              (*itEnd)=end;
              (*itBeg)=beg;
              ++itBeg;
              --itEnd;
              if(d==1)
                  break;
          }
          if (!d)
              flipCharacter(*itEnd, *itEnd);
      }

      void writeReverseComplementString(string &name, string &content, fstream &file)
      {
          file << ">" << name << "\n";
          flipNswap(content);
          file << content ;
      }

      void generateSeqPos(vector<seqData> &vecSeqInfo) {
          seqData s;
          uint64_t i=0,j=0;
          string line;
          clearFileFlag();
          if (numSeqFiles >= 1){
              while(getline(file1, line).good() ){
                  if(line[0] == '>'){
                      if(!strName.empty()) {
                          s.start=CHARS2BITS(j);
                          s.end=CHARS2BITS(i-1);
                          s.seq.assign(strtok(const_cast<char *>(strName.c_str())," \t\n"));
                          s.seq += "_1";
                          vecSeqInfo.push_back(s);
                          s.seq.clear();
                          i+=RANDOM_SEQ_SIZE;
                          j=i;
                          strName.clear();
                      }
                      if(!line.empty())
                          strName=line.substr(1);
                  } else if( !strName.empty() ) {
                      i+=line.length();
                  }
              }
              if (numSeqFiles == 2){
                  while(getline(file2, line).good() ){
                      if(line[0] == '>'){
                          if(!strName.empty()) {
                              s.start=CHARS2BITS(j);
                              s.end=CHARS2BITS(i-1);
                              s.seq.assign(strtok(const_cast<char *>(strName.c_str())," \t\n"));
                              s.seq += "_2";
                              vecSeqInfo.push_back(s);
                              s.seq.clear();
                              i+=RANDOM_SEQ_SIZE;
                              j=i;
                              strName.clear();
                          }
                          if(!line.empty())
                              strName=line.substr(1);
                      } else if( !strName.empty() ) {
                          i+=line.length();
                      }
                  }
              }
          }
          if( !strName.empty() ) {
              i+=line.length();
              s.start=CHARS2BITS(j);
              s.end=CHARS2BITS(i-1);
              s.seq.assign(strtok(const_cast<char *>(strName.c_str())," \t\n"));
              if (numSeqFiles == 1)
                  s.seq += "_1";
              if (numSeqFiles == 2)
                  s.seq += "_2";
              vecSeqInfo.push_back(s);
              s.seq.clear();
              strName.clear();
          }
      }

      void generateRevComplement() {
          string line,content;
          fstream revFile;

          char buffer[256];
          memset(buffer,0,256);
          sprintf(buffer, "%s/revComp", commonData::nucmer_path);
          revFile.open(buffer, ios::out);
          if (!revFile.is_open())
          {
              cout << "ERROR: unable to open temporary reverse complement file" << endl;
              exit( EXIT_FAILURE );
          }

          clearFileFlag();

          if (numSeqFiles >= 1){
              while(getline( file1, line ).good()){
                  size += (line.length()+1);
                  if(line[0] == '>'){
                      if(!strName.empty()) {
                          numSequences++;
                          size += RANDOM_SEQ_SIZE;
                          writeReverseComplementString(strName, content, revFile);
                          content.clear();
                          strName.clear();
                      }
                      if(!line.empty())
                          strName=line.substr(1);
                  } else if( !strName.empty() ) {
                      content += "\n";
                      content += line;
                  }
              }
              // Catches the last sequence line
              if( !strName.empty() ) {
                  size += (line.length()+1);
                  content += "\n";
                  content += line;
                  numSequences++;
                  writeReverseComplementString(strName, content, revFile);
                  content.clear();
                  strName.clear();
              }
              if (numSeqFiles == 2){
                  while(getline( file2, line ).good()){
                      size += (line.length()+1);
                      if(line[0] == '>'){
                          if(!strName.empty()) {
                              numSequences++;
                              size += RANDOM_SEQ_SIZE;
                              writeReverseComplementString(strName, content, revFile);
                              content.clear();
                              strName.clear();
                          }
                          if(!line.empty())
                              strName=line.substr(1);
                      } else if( !strName.empty() ) {
                          content += "\n";
                          content += line;
                      }
                  }
                  // Catches the last sequence line
                  if( !strName.empty() ) {
                      size += (line.length()+1);
                      content += "\n";
                      content += line;
                      numSequences++;
                      writeReverseComplementString(strName, content, revFile);
                      content.clear();
                      strName.clear();
                  }
              }
          }
          revFile.close();
      }
};

class outFileReadInfo {
    fstream outFile;

public:

    void setFile(string &filename)
    {
        outFile.open(filename, ios::out);
        if(!outFile.is_open()) {
            cout << "ERROR: unable to open " << filename << " file" << endl;
            exit( EXIT_FAILURE );
        }
    }


    void writeString(string &content, fstream &file)
    {
        file << content << "\n";
    }

    void removeSeq(vector<seqData> &vecSeqInfo, vector<string> &filenames) {
        seqData s;
        uint64_t lineNum = 0;
        string line;
        char buffer[256];
        memset(buffer,0,256);
        sprintf(buffer, "%s/IR", commonData::nucmer_path);
        fstream palFile;
        palFile.open(buffer, ios::in);

        /* Write inverted repeats */
        while (getline(palFile, line).good()) {
            writeString(line, outFile);
        }

        if (filenames.size() >= 1) {
            fstream file1;
            file1.open(filenames[0], ios::in);
            while (getline(file1, line).good()) {
                if (line[0] == '>') {
                    if (vecSeqInfo[lineNum].keep) {
                        line = strtok(const_cast<char *>(line.c_str()), " \t\n");
                        line += "_1";
                        writeString(line, outFile);
                    }
                } else {
                    if (vecSeqInfo[lineNum].keep)
                        writeString(line, outFile);
                    ++lineNum;
                }
            }
            file1.close();
            if (filenames.size() == 2) {
                fstream file2;
                file2.open(filenames[0], ios::in);
                while (getline(file2, line).good()) {
                    if (line[0] == '>') {
                        if (vecSeqInfo[lineNum].keep) {
                            line = strtok(const_cast<char *>(line.c_str()), " \t\n");
                            line += "_2";
                            writeString(line, outFile);
                        }
                    } else {
                        if (vecSeqInfo[lineNum].keep)
                            writeString(line, outFile);
                        ++lineNum;
                    }
                }
                file2.close();
            }
        }
    }

    void closeFile() {
        outFile.close();
    }
};

class SeqPos {
public:
    uint64_t lR;
    uint64_t rR;
    uint64_t lQ;
    uint64_t rQ;
    SeqPos() {
    }

    SeqPos(uint64_t lr, uint64_t rr, uint64_t lq, uint64_t rq)
    {
        lR=lr;
        rR=rr;
        lQ=lq;
        rQ=rq;
    }
};

class MemExt {
  public:
    uint64_t lR;
    uint64_t rR;
    uint64_t lQ;
    uint64_t rQ;
    MemExt() {
    }

    MemExt(uint64_t lr, uint64_t rr, uint64_t lq, uint64_t rq)
    {
         lR=lr;
         rR=rr;
         lQ=lq;
         rQ=rq;
    }

    bool operator () (const MemExt &obj1, const MemExt &obj2)
    {
      if (obj1.lQ<obj2.lQ)
          return true;
      else if (obj1.lQ>obj2.lQ)
          return false;
      else{
         if (obj1.lR<obj2.lR)
             return true;
         else
             return false;
      }
    }
};


class tmpFilesInfo {
    fstream *TmpFiles;
    fstream palFile;
    vector <MemExt> MemExtVec;
    uint64_t numMemsInFile;

    bool checkMEMExt(uint64_t &lr, uint64_t &rr, uint64_t &lq, uint64_t &rq, seqFileReadInfo &QueryFile, seqFileReadInfo &RefFile) {
      if ((!lq && QueryFile.getCurrPos()) || rq == CHARS2BITS(QueryFile.totalBases-1)) {
         return true;
      }else if((!lr && RefFile.getCurrPos()) || rr == CHARS2BITS(RefFile.totalBases-1)) {
         return true;
      }
      return false;
    }

    void writeToFile(uint64_t lQ, uint64_t rQ, uint64_t lR, uint64_t rR) {
        MemExt m;
        m.lQ=lQ;
        m.lR=lR;
        m.rQ=rQ;
        m.rR=rR;
        TmpFiles[m.lQ/numMemsInFile].write((char *)&m, sizeof(MemExt));
    }

    void writeToVector(uint64_t lQ, uint64_t rQ, uint64_t lR, uint64_t rR) {
        MemExt m;
        m.lQ=lQ;
        m.lR=lR;
        m.rQ=rQ;
        m.rR=rR;
        MemExtVec.emplace_back(m);
    }


  public:

    tmpFilesInfo(int numFiles) {
        TmpFiles = new fstream[numFiles];
    }

    ~tmpFilesInfo() {
        delete [] TmpFiles;
    }

    void setNumMemsInFile(uint64_t size, uint64_t &numSequences) {
        numMemsInFile = ((2*(size*commonData::d+numSequences*RANDOM_SEQ_SIZE+commonData::d))/NUM_TMP_FILES);
    }

    static bool compare_reference (const MemExt &obj1, const MemExt &obj2)
    {
      if (obj1.lR<obj2.lR)
         return true;
      else if (obj1.lR>obj2.lR)
         return false;
      else{
         if (obj1.lQ<obj2.lQ)
             return true;
         else
             return false;
      }

    }

    static bool sortReverse (const MemExt &obj1, const MemExt &obj2)
    {
        if (obj1.lQ + obj1.lR > obj2.lQ + obj2.lR)
            return true;
        else if (obj1.lQ + obj1.lR < obj2.lQ + obj2.lR)
            return false;
        else{
            if (obj1.rQ + obj1.rR > obj2.rQ + obj2.rR)
                return true;
            else
                return false;
        }
    }

    static bool myUnique(const MemExt &obj1, const MemExt &obj2)
    {
      if((obj1.lQ==obj2.lQ) && (obj1.rQ==obj2.rQ) && (obj1.rR==obj2.rR) && (obj1.lR==obj2.lR))
          return true;
      else
          return false;
    }

    static bool uniqueIR(const posData &obj1, const posData &obj2)
    {
        if((obj1.L1Bound==obj2.L2Bound) && (obj1.R1Bound==obj2.R1Bound) && (obj1.L2Bound==obj2.L2Bound) && (obj1.R2Bound==obj2.R2Bound))
            return true;
        else
            return false;
    }

    static bool myUniqueQue(const MemExt &obj1, const MemExt &obj2)
    {
        if((obj1.lQ==obj2.lQ) && (obj1.rQ==obj2.rQ))
            return true;
        else
            return false;
    }

    void openFiles(ios_base::openmode mode, int numFiles) {
        char buffer[256];
        memset(buffer,0,256);
        static int flag=0;

        sprintf(buffer, "%s", commonData::nucmer_path);
        if (!flag) {
            if(mkdir(buffer, S_IRWXU|S_IRGRP|S_IXGRP))
            {
                cout << "ERROR: unable to open temporary directory" << endl;
                exit( EXIT_FAILURE );
            }
            flag=1;
        }

        /* Last two files hold the sequence/pos mapping
         * for reference and query file respectively
         */
        for (int32_t i=0;i<numFiles;i++) {
            /* Temporary file to hold the mems */
            sprintf(buffer, "%s/%d", commonData::nucmer_path, i);
            TmpFiles[i].open(buffer, mode);
            if (!TmpFiles[i].is_open())
            {
                cout << "ERROR: unable to open temporary file" << endl;
                exit( EXIT_FAILURE );
            }
        }
    }

    void closeFiles(int numFiles) {
        for (int32_t i=0;i<numFiles;i++){
           TmpFiles[i].close();
        }
    }

    void setIRFile() {
        char buffer[256];
        memset(buffer,0,256);
        sprintf(buffer, "%s/IR", commonData::nucmer_path);

        palFile.open(buffer, ios::out);
        if(!palFile.is_open()) {
            cout << "ERROR: unable to open "<< buffer << " file" << endl;
            exit( EXIT_FAILURE );
        }
    }

    void closeIRFile() {
        palFile.close();
    }

    fstream& getMapFile(int fIndex) {
        return TmpFiles[fIndex];
    }

    bool writeMemInTmpFiles(uint64_t &lRef, uint64_t &rRef, uint64_t &lQue, uint64_t &rQue, seqFileReadInfo &QueryFile, seqFileReadInfo &RefFile) {
       MemExt m;
       uint64_t currPosQ = CHARS2BITS(QueryFile.getCurrPos());
       uint64_t currPosR = CHARS2BITS(RefFile.getCurrPos());
       if (rRef-lRef+2 >= static_cast<uint64_t>(commonData::minMemLen)) {
           if (!(commonData::d==1 && commonData::numThreads==1) && checkMEMExt(lRef, rRef, lQue, rQue, QueryFile, RefFile)) {
               #pragma omp critical(writeVector)
               writeToVector(currPosQ+lQue, currPosQ+rQue, currPosR+lRef,  currPosR+rRef);
           }else {
               #pragma omp critical(writeFile)
               writeToFile(currPosQ+lQue, currPosQ+rQue, currPosR+lRef,  currPosR+rRef);
           }
           return true;
       }else
           return false;
    }

    void printQueryHeader(vector<seqData>::iterator &itQ)
    {
        if (commonData::lenInHeader) {
            cout << "> " << (*itQ).seq << " Reverse" << " Len = " << ((*itQ).end-(*itQ).start+2)/2 << "\n";
        }else{
            cout << "> " << (*itQ).seq << " Reverse" << "\n";
        }
    }

    void printMemOnTerminal(vector<seqData> &refSeqInfo, vector<seqData> &querySeqInfo, MemExt &m) {
        uint64_t &lRef = m.lR;
        uint64_t &rRef = m.rR;
        uint64_t &lQue = m.lQ;
        uint64_t &rQue = m.rQ;
        static int flag=0;
        vector<seqData>::iterator itR;
        static vector<seqData>::iterator itQ=querySeqInfo.begin();
        seqData s;

        /* print remianing query sequences - if any */
        if (!lRef && !rRef && !lQue && !rQue) {
            /* No matches found - simply return */
            if (!flag){
                printQueryHeader(itQ);
            }
            while(itQ != --querySeqInfo.end()){
                ++itQ;
                printQueryHeader(itQ);
            }
            itQ=querySeqInfo.begin();
            flag=0;
            return;
        }


        s.start=lRef;
        s.end=rRef;

        if (rRef-lRef+2 < static_cast<uint64_t>(commonData::minMemLen))
            return;


        /* Process relative position for Reference sequence */
        itR = lower_bound(refSeqInfo.begin(), refSeqInfo.end(), s, seqData());

        if ((*itR).start <= lRef && (*itR).end >= rRef){
            // MEM within acutal sequence
            // s------e--s------e--s------e
            // s--|--|e
            lRef?lRef-=((*itR).start):lRef;
            rRef-=((*itR).start);
        }else if ((*itR).start > lRef && (*itR).end >= rRef) {
            if ((*itR).start > rRef) //mem within random character
                 return;
            // s------e--s------e--s------e
            // s------e-|s--|---e
            lQue+=((*itR).start-lRef);
            lRef=0;
            rRef-=((*itR).start);
        }else if ((*itR).start > lRef && (*itR).end < rRef) {
            // s------e--s------e--s------e
            // s------e-|s------e-|s------e
            lQue+=((*itR).start-lRef);
            lRef=0;
            rQue-=(rRef-(*itR).end);
            rRef=((*itR).end-(*itR).start);
        }else if ((*itR).start <= lRef && (*itR).end < rRef) {
            // s------e--s------e--s------e
            // s------e--s-----|e-|s------e
            rQue-=(rRef-(*itR).end);
            lRef?lRef-=((*itR).start):lRef;
            rRef=((*itR).end-(*itR).start);
        }else //mem within random character
            return;

        if (rRef-lRef+2 < static_cast<uint64_t>(commonData::minMemLen))
            return;

        /* Print first Query sequence */
        if (!flag){
            printQueryHeader(itQ);
            flag=1;
        }
        /* Process relative position for Query sequence */
        while(lQue >= (*itQ).end){
            ++itQ;
            printQueryHeader(itQ);
        }
        if ((*itQ).start <= lQue && (*itQ).end >= rQue){
            // MEM within acutal sequence
            // s------e--s------e--s------e
            // s--|--|e
            lQue?lQue-=((*itQ).start):lQue;
            rQue-=((*itQ).start);
        }else if ((*itQ).start > lQue && (*itQ).end >= rQue) {
            if ((*itQ).start > rQue) //mem within random character
                 return;
            // s------e--s------e--s------e
            // s------e-|s--|---e
            lRef+=((*itQ).start-lQue);
            lQue=0;
            rQue-=((*itQ).start);
        }else if ((*itQ).start > lQue && (*itQ).end < rQue) {
            // s------e--s------e--s------e
            // s------e-|s------e-|s------e
            lRef+=((*itQ).start-lQue);
            lQue=0;
            rRef-=(rQue-(*itQ).end);
            rQue=((*itQ).end-(*itQ).start);
        }else if ((*itQ).start <= lQue && (*itQ).end < rQue) {
            // s------e--s------e--s------e
            // s------e--s-----|e-|s------e
            rRef-=(rQue-(*itQ).end);
            lQue?lQue-=((*itQ).start):lQue;
            rQue=((*itQ).end-(*itQ).start);
        }else //mem within random character
            return;


        if (rRef-lRef+2 >= static_cast<uint64_t>(commonData::minMemLen)){
           if (refSeqInfo.size() == 1 && !commonData::fourColOutput) {
               if (commonData::relQueryPos)
                   cout << " " << setw(15) << ((lRef+2)/2) <<  setw(15) << ((*itQ).end-(*itQ).start-lQue+2)/2 << setw(15) << ((rRef-lRef+2)/2) << "\n";
               else
                   cout << " " << setw(15) << ((lRef+2)/2) <<  setw(15) << ((lQue+2)/2) << setw(15) << ((rRef-lRef+2)/2) << "\n";
           }else{
               if (commonData::relQueryPos) {
                   cout << " " << setw(30) << std::left <<(*itR).seq << setw(15) << ((lRef+2)/2) <<  setw(15) << ((*itQ).end-(*itQ).start-lQue+2)/2 << setw(15) << ((rRef-lRef+2)/2) << "\n";
               }else{
                   cout << " " << setw(30) << std::left <<(*itR).seq << setw(15) << ((lRef+2)/2) <<  setw(15) << ((lQue+2)/2) << setw(15) << ((rRef-lRef+2)/2) << "\n";
               }
           }
        }
    }

    void mergeMemExtVector () {
        int flag=0;
        MemExt m;
        //if (commonData::d==1 && commonData::numThreads==1)
        //    return;

        if (MemExtVec.size() > 1) {
            sort(MemExtVec.begin(), MemExtVec.end(), MemExt());
            //Remove duplicates, this is happens due to removal of all checks.
            vector<MemExt>::iterator last=unique(MemExtVec.begin(), MemExtVec.end(), myUnique);
            do {
                flag=0;
                for (vector<MemExt>::iterator it=MemExtVec.begin(); it != (last-1); ++it) {
                    if (!(*it).lQ && !(*it).rQ && !(*it).lR && !(*it).rR )
                        continue;

                    vector<MemExt>::iterator dup = it;
                    ++dup;
                     for (; dup != last; ++dup) {
                        if (!(*dup).lQ && !(*dup).rQ && !(*dup).lR && !(*dup).rR )
                            continue;
                        if((*dup).lQ + static_cast<uint64_t>(commonData::minMemLen-2)-2 > (*it).rQ)
                            break;
                        if((*dup).lQ + static_cast<uint64_t>(commonData::minMemLen-2)-2 == (*it).rQ) {
                            if((*dup).lR + static_cast<uint64_t>(commonData::minMemLen-2)-2 == (*it).rR) {
                                flag=1;
                                (*it).rQ=(*dup).rQ;
                                (*it).rR=(*dup).rR;
                                (*dup).rQ=0;
                                (*dup).rR=0;
                                (*dup).lQ=0;
                                (*dup).lR=0;
                                break;
                            }
                        }
                    }
                    if (flag)
                        break;
                }
            } while (flag);

            last=unique(MemExtVec.begin(), last, myUnique);

            sort(MemExtVec.begin(), last, compare_reference);
            do {
                flag=0;
                for (vector<MemExt>::iterator it=MemExtVec.begin(); it != (last-1); ++it) {
                    if (!(*it).lQ && !(*it).rQ && !(*it).lR && !(*it).rR )
                        continue;
                    vector<MemExt>::iterator dup = it;
                    ++dup;
                    for (; dup != last; ++dup) {
                        if (!(*dup).lQ && !(*dup).rQ && !(*dup).lR && !(*dup).rR )
                            continue;
                        if((*dup).lR + static_cast<uint64_t>(commonData::minMemLen-2)-2 > (*it).rR)
                            break;
                        if((*dup).lR + static_cast<uint64_t>(commonData::minMemLen-2)-2 == (*it).rR) {
                            if((*dup).lQ + static_cast<uint64_t>(commonData::minMemLen-2)-2 == (*it).rQ) {
                                flag=1;
                                (*it).rQ=(*dup).rQ;
                                (*it).rR=(*dup).rR;
                                (*dup).rQ=0;
                                (*dup).rR=0;
                                (*dup).lQ=0;
                                (*dup).lR=0;
                                break;
                            }
                        }
                    }
                    if (flag)
                        break;
                }
            } while (flag);

            for (vector<MemExt>::iterator it=MemExtVec.begin(); it != last; ++it) {
                if ((*it).lQ || (*it).rQ || (*it).lR || (*it).rR ) {
                    writeToFile((*it).lQ, (*it).rQ, (*it).lR, (*it).rR);
                }
            }
        }
        MemExtVec.clear();
    }

    void extBin(seqFileReadInfo &RefFile, vector<uint64_t> &Bins, int &binLength, uint64_t &extLength, uint64_t &Bound) {
        int readBin=0, bin=0;
        uint64_t offset, totalOffset=0;

        Bins.resize(binLength);
        while (bin != binLength) {
            offset = (Bound) % DATATYPE_WIDTH;
            if (totalOffset % DATATYPE_WIDTH) {
                offset = DATATYPE_WIDTH - offset;
                Bins[bin] |= ((RefFile.binReads[(Bound) / DATATYPE_WIDTH + readBin] & global_mask_left[(totalOffset % DATATYPE_WIDTH)/2 - 1]) >> offset);
                ++bin;
            } else {
                Bins[bin] = RefFile.binReads[(Bound) / DATATYPE_WIDTH + readBin];
                Bins[bin] <<= offset;
                ++readBin;
                if (!offset)
                    ++bin;
            }
            totalOffset += offset;
        }
        if (extLength % DATATYPE_WIDTH)
            Bins[binLength-1] = Bins[binLength-1] & global_mask_left[(extLength % DATATYPE_WIDTH)/2 - 1];
    }

    void convertToNucl(uint64_t bin, string &sequence) {
        if (bin == 0) {
            sequence += 'A';
        } else if (bin == 1) {
            sequence += 'C';
        } else if (bin == 2) {
            sequence += 'G';
        } else if (bin == 3) {
            sequence += 'T';
        }
    }

    int hammingDistance(uint64_t &bin1, uint64_t &bin2)
    {
        uint64_t x = bin1 ^ bin2;
        int setBits = 0;

        while (x > 0) {
            setBits += x & 1;
            x >>= 1;
        }

        return setBits;
    }

    void writeInvertedRepeats(seqFileReadInfo &RefFile, posData &posDataInfo) {

        uint64_t currBin;
        uint32_t offset;

        string sequence;
        for (uint64_t pos = ((posDataInfo.L1Bound) / 2); pos != ((posDataInfo.R1Bound) / 2); ++pos) {
            currBin = RefFile.binReads[(pos * 2) / DATATYPE_WIDTH];
            offset = (pos * 2) % DATATYPE_WIDTH;
            currBin &= global_mask_right[(DATATYPE_WIDTH - offset) / 2 - 1];
            currBin >>= ((DATATYPE_WIDTH - 2) - offset);
            convertToNucl(currBin, sequence);
        }
        for (uint64_t pos = ((posDataInfo.L2Bound) / 2); pos != ((posDataInfo.R2Bound) / 2); ++pos) {
            currBin = RefFile.binReads[(pos * 2) / DATATYPE_WIDTH];
            offset = (pos * 2) % DATATYPE_WIDTH;
            currBin &= global_mask_right[(DATATYPE_WIDTH - offset) / 2 - 1];
            currBin >>= ((DATATYPE_WIDTH - 2) - offset);
            convertToNucl(currBin, sequence);
        }
        palFile << posDataInfo.seq << "_LCoord_" << ((posDataInfo.R1Bound - posDataInfo.L1Bound) / 2 - (posDataInfo.palLength / 2) + 1) << "_RCoord_" << ((posDataInfo.R1Bound - posDataInfo.L1Bound) / 2) << "\n";
        palFile << sequence << "\n";
    }

    void getInvertedRepeats(seqFileReadInfo &RefFile, vector<seqData> &vecSeqInfo) {
        int numFiles=0;
        vector<MemExt> MemExtVec;
        MemExt m;
        vector<uint64_t> l1Bins, l2Bins, r1Bins, r2Bins;
        mapObject QueryNpos, RefNpos1, RefNpos2;
        uint64_t Bound, currL1Bound, currR1Bound, currL2Bound, currR2Bound, ext=2, extLengthL, extLengthR, currExtLengthL, currExtLengthR;
        uint64_t duplR = 0, duprR = 0, duplQ = 0, duprQ = 0;
        int binLength, hDL, hDR;
        seqData s;
        vector<seqData>::iterator seqit;
        string currHeader;
        posData p;
        int flag = 1;
        char buffer[256];
        memset(buffer,0,256);

        numFiles=NUM_TMP_FILES;

        //sprintf(buffer, "%s/%d", commonData::nucmer_path, numFiles);
        //remove(buffer);

        //sprintf(buffer, "%s/%d", commonData::nucmer_path, numFiles+1);
        //remove(buffer);

        openFiles(ios::in|ios::binary, numFiles);
        setIRFile();

        for (int32_t i=0;i<numFiles;i++) {
            sprintf(buffer, "%s/%d", commonData::nucmer_path, i);
            if (i == NUM_TMP_FILES) {
                /* Output any unsued query sequence */
                m.lR = m.lQ = m.rR = m.rQ = 0;
                /* Redirect output to reverse complement file */
                std::cout.rdbuf(TmpFiles[numFiles + 1].rdbuf());
            }
            while (!TmpFiles[i].read((char *) &m, sizeof(MemExt)).eof()) {
                MemExtVec.emplace_back(m);
            }
            TmpFiles[i].close();
        }

        sort(MemExtVec.begin(), MemExtVec.end(), myUniqueQue);

        for (vector<MemExt>::iterator it=MemExtVec.begin();it!=MemExtVec.end();++it) {

            if (flag) {
                RefFile.getKmerLeftnRightBoundForNs((*it).lR, RefNpos1);
                currExtLengthL = (*it).lR - RefNpos1.left;
                currExtLengthR = RefNpos1.right - (*it).rR;
                currL1Bound = RefNpos1.left;
                currR1Bound = (*it).rR + ext;
                currL2Bound = (*it).rR + ext;
                currR2Bound = RefNpos1.right + ext;
            }

            if (duprQ == (*it).rQ && duplQ == (*it).lQ) {
                RefFile.getKmerLeftnRightBoundForNs(duplR, RefNpos1);
                RefFile.getKmerLeftnRightBoundForNs((*it).lR, RefNpos2);
                // Left extension
                if (max(duplR - RefNpos1.left, (*it).lR - RefNpos2.left) > currExtLengthL) {
                    currExtLengthL = max(duplR - RefNpos1.left, (*it).lR - RefNpos2.left);

                    extLengthL = min(duplR - RefNpos1.left, (*it).lR - RefNpos2.left);
                    binLength = min((duplR - RefNpos1.left) / DATATYPE_WIDTH + 1,
                                    ((*it).lR - RefNpos2.left) / DATATYPE_WIDTH + 1);
                    Bound = duplR - extLengthL;
                    extBin(RefFile, l1Bins, binLength, extLengthL, Bound);
                    Bound = (*it).lR - extLengthL;
                    extBin(RefFile, l2Bins, binLength, extLengthL, Bound);

                    hDL = 0;
                    for (int bin = 0; bin != binLength; ++bin) {
                        if (l1Bins[bin] != l2Bins[bin])
                            hDL += hammingDistance(l1Bins[bin], l2Bins[bin]);
                    }

                    if (hDL <= HAMMING_DISTANCE_LIM) {
                        if (duplR - RefNpos1.left >= (*it).lR - RefNpos2.left) {
                            currL1Bound = RefNpos1.left;
                            currR1Bound = duprR + ext;
                        } else {
                            currL1Bound = RefNpos2.left;
                            currR1Bound = (*it).rR + ext;
                        }
                    }
                }

                // Right extension
                if (max(RefNpos1.right - duprR, RefNpos2.right - (*it).rR) > currExtLengthR) {
                    currExtLengthR = max(RefNpos1.right - duprR, RefNpos2.right - (*it).rR);

                    extLengthR = min(RefNpos1.right - duprR, RefNpos2.right - (*it).rR);
                    binLength = min((RefNpos1.right - duprR) / (DATATYPE_WIDTH + 1) + 1,
                                    (RefNpos2.right - (*it).rR) / (DATATYPE_WIDTH + 1) + 1);
                    Bound = duprR + ext;
                    extBin(RefFile, r1Bins, binLength, extLengthR, Bound);
                    Bound = (*it).rR + ext;
                    extBin(RefFile, r2Bins, binLength, extLengthR, Bound);

                    hDR = 0;
                    for (int bin = 0; bin != binLength; ++bin) {
                        if (r1Bins[bin] != r2Bins[bin])
                            hDR += hammingDistance(r1Bins[bin], r2Bins[bin]);
                    }

                    if (hDR <= HAMMING_DISTANCE_LIM) {
                        if (RefNpos1.right - duprR >= RefNpos2.right - (*it).rR) {
                            currL2Bound = duprR + ext;
                            currR2Bound = RefNpos1.right + ext;
                        } else {
                            currL2Bound = (*it).rR + ext;
                            currR2Bound = RefNpos2.right + ext;
                        }
                    }
                }

                if (flag) {
                    s.start=duprR;
                    s.end=duplR;
                    seqit = lower_bound(vecSeqInfo.begin(), vecSeqInfo.end(), s, seqData());
                    currHeader += "_from_";
                    currHeader += (*seqit).seq;
                    (*seqit).keep = 0;
                }

                s.start=(*it).rR;
                s.end=(*it).lR;
                seqit = lower_bound(vecSeqInfo.begin(), vecSeqInfo.end(), s, seqData());
                currHeader += "_from_";
                currHeader += (*seqit).seq;
                (*seqit).keep = 0;

                flag = 0;
            } else {

                if (!flag) {
                    p.L1Bound=currL1Bound;
                    p.R1Bound=currR1Bound;
                    p.L2Bound=currL2Bound;
                    p.R2Bound=currR2Bound;
                    p.palLength=((*it).rR - (*it).lR + ext);
                    p.seq=currHeader;
                    writeInvertedRepeats(RefFile, p);
                    flag = 1;
                }

                duplR = (*it).lR;
                duplQ = (*it).lQ;
                duprR = (*it).rR;
                duprQ = (*it).rQ;
                currHeader = ">InvertedRepeat";
                s.start=duprQ;
                s.end=duplQ;
                seqit = lower_bound(vecSeqInfo.begin(), vecSeqInfo.end(), s, seqData());
                currHeader += "_of_";
                currHeader += (*seqit).seq;
            }
        }
        closeIRFile();
    }

    void removeTmp() {
        char buffer[256];
        memset(buffer,0,256);

        sprintf(buffer, "%s/revComp", commonData::nucmer_path);
        remove(buffer);

        sprintf(buffer, "%s/IR", commonData::nucmer_path);
        remove(buffer);

        sprintf(buffer, "%s", commonData::nucmer_path);
        remove(buffer);

    }
};
