/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include <iostream>
#include <map>

#include "maptiles.h"

using namespace std;


Point<3> convertToXYZ(LUVAPixel pixel) {
    return Point<3>( pixel.l, pixel.u, pixel.v );
}

MosaicCanvas* mapTiles(SourceImage const& theSource,
                       vector<TileImage>& theTiles)
{
    vector<Point<3>> converted_tiles;
    map<Point<3>, TileImage*> Converted;
    for (auto& pt : theTiles) {
        LUVAPixel pixel = pt.getAverageColor();
        Point<3> point = convertToXYZ(pixel);
        converted_tiles.push_back(point);
        Converted[point] = &pt;
    }
    KDTree<3> AvgColourTree(converted_tiles);
    MosaicCanvas *ReturnCanvas = new MosaicCanvas(theSource.getRows(), theSource.getColumns());
    for (int y = 0; y < ReturnCanvas->getColumns(); y++) {
        for (int x = 0; x < ReturnCanvas->getRows(); x++) {
            Point<3> cur = convertToXYZ(theSource.getRegionColor(x, y));
            Point<3> match = AvgColourTree.findNearestNeighbor(cur);
            ReturnCanvas->setTile(x, y, Converted[match]);
        }
    }
    return ReturnCanvas;
}

