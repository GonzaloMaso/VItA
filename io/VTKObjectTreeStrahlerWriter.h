#ifndef VTKOBJECTTREESTRAHLERWRITER_H_
#define VTKOBJECTTREESTRAHLERLWRITER_H_

#include "VTKObjectTreeWriter.h"

class VTKObjectTreeStrahlerWriter: public VTKObjectTreeWriter {
public:
	VTKObjectTreeStrahlerWriter();
	virtual void write(string filename, AbstractObjectCCOTree *tree);
	virtual ~VTKObjectTreeStrahlerWriter();
};

#endif /* VTKOBJECTTREESTRAHLERWRITER_H_ */