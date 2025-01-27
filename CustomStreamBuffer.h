#pragma once
#include <streambuf>
#include <iostream>
#include <sstream>
#include <FECore/log.h>
#include <FECore/FEModel.h>

class CustomStreamBuffer : public std::streambuf 
{
public:
   CustomStreamBuffer() {};
   CustomStreamBuffer(std::streambuf* buf, FEModel *fe_model_) : originalBuffer(buf), fe_model(fe_model_) {};

   FEModel* GetFEModel() const { return fe_model; };

   std::streambuf* originalBuffer = nullptr;
   FEModel* fe_model = nullptr;
   
protected:
   virtual int overflow(int c) override {
      if (c != EOF) {
         char ch = static_cast<char>(c);
         buffer << ch;
         if (ch == '\n') {
               feLog("%s", buffer.str().c_str());
               buffer.str("");
               buffer.clear();
         }
      }
      return c;
   }

private:
   std::ostringstream buffer;
};