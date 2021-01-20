// Equivalent code to mdat_conv.py, but implemented as a ROOT macro to avoid having to build python3 ROOT libraries



// -------------------------------------------------------------------//
// ------------------------- Global variables ------------------------//
// -------------------------------------------------------------------//

// Masks for extracting amp, pos etc. from 48 bit event block
ULong64_t mask_eventID = 0b100000000000000000000000000000000000000000000000;
ULong64_t mask_amp =     0b011111111000000000000000000000000000000000000000;
ULong64_t mask_ypos =    0b000000000111111111100000000000000000000000000000;
ULong64_t mask_xpos =    0b000000000000000000011111111110000000000000000000;
ULong64_t mask_time =    0b000000000000000000000000000001111111111111111111;


// Header parameters
struct Header{
  uint16_t bufferlength;
  uint16_t buffertype;
  uint16_t headerlength;
  uint16_t buffernumber;
  uint16_t runID;
  uint8_t  mcpdID;          // 1 for segment 1, 2 for segment 2
  uint8_t  status;          //
  uint64_t headerTS;        // 48 bit timestamp
  uint64_t param0;          // 48 bit parameter - unused
  uint64_t param1;          // 48 bit parameter - unused
  uint64_t param2;          // 48 bit parameter - unused
  uint64_t param3;          // 48 bit parameter - unused
} header;

// Individual event parameters
struct Event{
  uint16_t xpos;           // wire number 
  uint16_t ypos;           // stripe number
  uint16_t amp;            // ToT in clock cycles (12.5 ns)
  uint64_t time;           // The full time stamp in clocks (12.5 ns)
  uint8_t eventID;         // 0 for real events, 1 for self triggers 
  uint32_t eventTS;        // The 19 bit time stamp within the buffer
} event;


// -------------------------------------------------------------------//
// --------------------- File reading functions ----------------------//
// -------------------------------------------------------------------//

// Byteswap a two-byte word
void ByteSwap16(uint16_t &word){
  word = (word)<<8 | (word)>>8;
}

// Read in a single two-byte word and byteswap (calling byteswap function)
void ReadWord(ifstream &infile, uint16_t &word){
  infile.read((char*) &word, 2);
  ByteSwap16(word);
}

// Read in a single byte
void ReadByte(ifstream &infile, uint8_t &byte){
  infile.read((char*) &byte, 1);
}

// Read a six-byte event, consisting of 3 words
// Used for event data and some header parameters
void ReadEntry(ifstream &infile, uint64_t &entry){
  uint16_t low, mid, high;
  ReadWord(infile, low);
  ReadWord(infile, mid);
  ReadWord(infile, high);
  entry = (uint64_t)low | (uint64_t)mid<<16 | (uint64_t)high<<32;
} 

// Read (and dispose of) 58 bytes of file header
void ReadHeader(ifstream &infile){ 
  char buffer[58];
  infile.read(buffer, 58);
}

// Read in an event buffer
int ReadBuffer(ifstream &infile){

  ReadWord(infile, header.bufferlength);  // Length of the buffer
  ReadWord(infile, header.buffertype);    // Should be 0x0002
  if(header.buffertype != 0x0002){
    return 1;
  }
  ReadWord(infile, header.headerlength);  // Should be constant
  ReadWord(infile, header.buffernumber);
  ReadWord(infile, header.runID); 
  ReadByte(infile, header.mcpdID);
  ReadByte(infile, header.status);
  ReadEntry(infile, header.headerTS);
  ReadEntry(infile, header.param0);
  ReadEntry(infile, header.param1);
  ReadEntry(infile, header.param2);
  ReadEntry(infile, header.param3);

  return 0;

}

// Read a single 48 bit event and split into the component parts
void ReadEvent(ifstream &infile){
  uint64_t rawevent;
  ReadEntry(infile, rawevent);
  event.eventID = (rawevent & mask_eventID) >> 47;
  event.amp = (rawevent & mask_amp) >> 39;
  event.ypos = (rawevent & mask_ypos) >> 29;
  event.xpos = (rawevent & mask_xpos) >> 19;
  event.eventTS = (rawevent & mask_time);
  event.time = event.eventTS + header.headerTS;
}

void ReadBufferEnd(ifstream &infile, int debug){
  uint16_t word;
  if ((debug & 4) > 0) cout << "--- Buffer padding ---" << endl;
  for (int i = 0; i < 4; i++){
    ReadWord(infile, word);
    if ((debug & 4) > 0) cout << hex << word << endl;
  }
  
}

// Print the current event buffer
void PrintBuffer(){
  cout << "----------------------------------------------------" << endl;
  cout << "Buffer number: " << header.buffernumber << endl;
  cout << "Buffer length: " << header.bufferlength << endl;
  cout << "Expected number of entries: " << (header.bufferlength -21)/3 << endl;
  cout << "Header length: " << header.headerlength << endl;
  cout << "Run ID: " << header.runID << endl;
  cout << "MCPD ID: " << int(header.mcpdID) << endl;
  cout << "Status: " << int(header.status) << endl;
  cout << "Header timestamp: " << header.headerTS << endl;
  cout << "Parameter 0: " << header.param0 << endl;
  cout << "Parameter 1: " << header.param1 << endl;
  cout << "Parameter 2: " << header.param2 << endl;
  cout << "Parameter 3: " << header.param3 << endl;
  cout << "----------------------------------------------------" << endl;
}

// Print the current event information
void PrintEvent(){
  cout << "----------------------------------------------------" << endl;
  cout << "EventID: " << int(event.eventID) << endl;
  cout << "xpos: " << event.xpos << endl;
  cout << "ypos: " << event.ypos << endl;
  cout << "amp: " << event.amp << endl;
  cout << "time stamp: " << event.eventTS << endl;
  cout << "absolute time: " << event.time << endl;
  cout << "----------------------------------------------------" << endl;
}

// -------------------------------------------------------------------//
// ------------------------------ Main -------------------------------//
// -------------------------------------------------------------------//


// debug 0 = off, 1 = buffer, 2 = events, 4 = post-buffer padding, 7 = all
void mdat_conv(TString filename, int debug=0){

  uint64_t total_events = 0;  // Total event counter
  uint64_t buffer_num = 0;    // Current buffer number
  uint64_t entry_num = 0;     // Current entry (event) number
  
  //--- Create the output ROOT file ---//

  // Set output filename
  TString outfilename = filename;
  outfilename.ReplaceAll(".mdat",".root");  
  
  TFile *outfile = new TFile(outfilename,"RECREATE");
  TTree *rawdata = new TTree("rawdata","Raw data converted from mdat to ROOT");

  // Set up branches - link to struct parameters
  rawdata->Branch("xpos", &event.xpos, "xpos/s");
  rawdata->Branch("ypos", &event.ypos, "ypos/s");
  rawdata->Branch("amp", &event.amp, "amp/s");
  rawdata->Branch("time", &event.time, "time/l");
  rawdata->Branch("eventID", &event.eventID, "eventID/b");
  rawdata->Branch("eventTS", &event.eventTS, "eventTS/i");
  rawdata->Branch("mcpdID", &header.mcpdID, "mcpdID/b");
  rawdata->Branch("status", &header.status, "status/b");
  rawdata->Branch("param0", &header.param0, "param0/l");
  rawdata->Branch("param1", &header.param1, "param1/l");
  rawdata->Branch("param2", &header.param2, "param2/l");
  rawdata->Branch("param3", &header.param3, "param3/l");  
  


  //--- Open the input mdat file ---//
  
  ifstream infile(filename, ios::in | ios::binary); 
  
  // Read past the file header, 58 bytes
  ReadHeader(infile);
  
  // Loop over all buffers - break loop when the wrong buffer header type is found
  
  while(true){
  
    if (ReadBuffer(infile) == 1) break;
    else buffer_num++;
    
    if ((debug & 1) > 0) PrintBuffer();
    
    int buffer_entries = (header.bufferlength - 21) / 3; // Expected number of entries in the current buffer
    
    for (int bentry = 0; bentry < buffer_entries; bentry++){
      
      ReadEvent(infile);
      if ((debug & 2) > 0) PrintEvent();
      entry_num++;
      rawdata->Fill();
      
      // Print info on status
      if (entry_num % 10000 == 0){
        cout << "Processing entry number: " << entry_num << "\r" << flush;
      }
    }
    
    // Read past the end-of-buffer words
    ReadBufferEnd(infile, debug);
  
  }
  
  cout << "---------------------------------------------------------" << endl;
  cout << "A total of " << entry_num << " events were read from " << buffer_num << " buffers" << endl; 
  cout << "---------------------------------------------------------" << endl;

  outfile->Write();

  // Close files
  infile.close();
  outfile->Close();
}
