/* SPDX-License-Identifier: Apache-2.0 */
/* Copyright 2020 Gonzalo Maso Talou */
/*
 * EXFTreeWriter.cpp
 *
 *  Created on: 27/07/2018
 *      Author: Gonzalo D. Maso Talou
 */

//#include "EXFTreeWriter.h"
//
//EXFTreeWriter::EXFTreeWriter(){
//}
//
//EXFTreeWriter::~EXFTreeWriter(){
//}
//
//void EXFTreeWriter::writeTree(string filename, AbstractStructuredCCOTree* tree)	{
//
//	vector< vector<double> > points; 		/** Triple of coordinates describing a point */
//	vector< vector<double> > directions; 	/** Triple of coordinates the point's direction */
//	vector< vector<int> > lines;			/** Tuple of points describing a line */
//	vector< vector<double> > radii;			/**	Single value describing the radius per vessel */
//
//	//	Computing topology
//	int idPoints = 0l;
//	vector<AbstractVascularElement *> elements = tree->getSegments();
//	for (std::vector<AbstractVascularElement *>::iterator it = elements.begin(); it != elements.end(); ++it) {
//		for (vector<SingleVessel *>::iterator it2 = (*it)->getVessels().begin(); it2 != (*it)->getVessels().end(); ++it2) {
//			SingleVessel *currentSegment = *it2;
//
//			int idProx = ++idPoints;
//			vector<double> pointProx;
//			pointProx.assign(currentSegment->xProx.p, currentSegment->xProx.p + 3);
//			points.push_back(pointProx);
//
//			int idDist = ++idPoints;
//			vector<double> pointDist;
//			pointDist.assign(currentSegment->xDist.p, currentSegment->xDist.p + 3);
//			points.push_back(pointDist);
//
//			point dirParent = {0.,0.,0.};
//			if(currentSegment->parent)
//				dirParent = ((SingleVessel*)currentSegment->parent)->xDist - ((SingleVessel*)currentSegment->parent)->xProx;
//			point dirVessel = currentSegment->xDist - currentSegment->xProx;
//			dirParent = (dirParent + dirVessel) * 0.5;
//			directions.push_back({ dirParent.p[0], dirParent.p[1], dirParent.p[2]}); 	// Proximal direction
//			directions.push_back({ dirVessel.p[0], dirVessel.p[1], dirVessel.p[2]});	// Distal direction
//
//			lines.push_back({idProx, idDist});
//		}
//	}
//
//
//	//	Reading tree information
//	for (std::vector<AbstractVascularElement *>::iterator it = elements.begin(); it != elements.end(); ++it) {
//		for (vector<SingleVessel *>::iterator it2 = (*it)->getVessels().begin(); it2 != (*it)->getVessels().end(); ++it2) {
//			SingleVessel *currentSegment = *it2;
//			vector<double> radius;
//
//			radius.push_back(currentSegment->radius);
//			radii.push_back(radius);
//		}
//	}
//
//	//	OpenCMISS structures
//	Context context("linear");
//	Region region = context.getDefaultRegion();
//	Fieldmodule fm = region.getFieldmodule();
//
//	fm.beginChange();
//
//	//	Declare fields to be defined on nodes and elements below
//	FieldFiniteElement coordinates = fm.createFieldFiniteElement(3);
//	coordinates.setName("coordinates");
//	coordinates.setManaged(true);
//	coordinates.setTypeCoordinate(true);
//	coordinates.setCoordinateSystemType(Field::COORDINATE_SYSTEM_TYPE_RECTANGULAR_CARTESIAN);
//	coordinates.setComponentName(1,"x");
//	coordinates.setComponentName(2,"y");
//	coordinates.setComponentName(3,"z");
//
//	Nodeset nodes = fm.findNodesetByFieldDomainType(Field::DOMAIN_TYPE_NODES);
//	Nodetemplate nodeTemplate = nodes.createNodetemplate();
//	nodeTemplate.defineField(coordinates);
//	nodeTemplate.setValueNumberOfVersions(coordinates, -1, Node::VALUE_LABEL_D_DS1, 1);
//
//	Fieldcache cache = fm.createFieldcache();
//
//	for (unsigned int n = 0; n < points.size(); ++n) {
//		Node node = nodes.createNode(n+1, nodeTemplate);
//		cache.setNode(node);
////		coordinates.assignReal(cache, 3, &(points[n][0]));
//		coordinates.setNodeParameters(cache,-1,Node::VALUE_LABEL_VALUE,1,3,&(points[n][0]));
//		coordinates.setNodeParameters(cache,-1,Node::VALUE_LABEL_D_DS1,1,3,&(directions[n][0]));
//	}
//
//	FieldFiniteElement radius = fm.createFieldFiniteElement(1);
//	radius.setName("radius");
//	radius.setManaged(true);
//	radius.setComponentName(1, "value");
//
//	Mesh mesh = fm.findMeshByDimension(1);
//	Elementtemplate elementTemplate = mesh.createElementtemplate();
//	elementTemplate.setElementShapeType(Element::SHAPE_TYPE_LINE);
//	Elementbasis cubicHermiteBasis = fm.createElementbasis(1, Elementbasis::FUNCTION_TYPE_CUBIC_HERMITE);
//	Elementfieldtemplate eftHermite = mesh.createElementfieldtemplate(cubicHermiteBasis);
//	elementTemplate.defineField(coordinates,-1,eftHermite);
////	Elementbasis linearBasis = fm.createElementbasis(1,Elementbasis::FUNCTION_TYPE_LINEAR_LAGRANGE);
////	Elementfieldtemplate eftLinear = mesh.createElementfieldtemplate(linearBasis);
////	elementTemplate.defineField(coordinates,-1,eftLinear);
//	Elementbasis constantBasis = fm.createElementbasis(1,Elementbasis::FUNCTION_TYPE_CONSTANT);
//	Elementfieldtemplate eftConstant = mesh.createElementfieldtemplate(constantBasis);
//	eftConstant.setParameterMappingMode(Elementfieldtemplate::PARAMETER_MAPPING_MODE_ELEMENT);
//	elementTemplate.defineField(radius,-1,eftConstant);
//
//	for (unsigned int e = 0; e < lines.size(); ++e) {
//		Element element = mesh.createElement(e + 1, elementTemplate);
//		element.setNodesByIdentifier(eftHermite,2,&(lines[e][0]));
//		cache.setElement(element);
//		radius.assignReal(cache, 1, &(radii[e][0]));
//	}
//
//	fm.endChange();
//
//	region.writeFile(filename);
//}
