//  Matrix.c
//  ModVoc
//
//  Created by Edward Costello on 23/08/2014.
//  Copyright (c) 2014 Edward Costello. All rights reserved.
//
//

#import "Matrix.h"
#import <string.h>

static bool staticTrue = true;
static bool staticFalse = false;

#pragma mark - Allocation / Deallocation -

OVERLOADED Matrix *Matrix_new(size_t rowCount , size_t columnCount)
{
    Matrix *self = calloc(1, sizeof(Matrix));
    self->rowCount = rowCount;
    self->columnCount = columnCount;
    self->rowStride = columnCount;
    self->columnStride = rowCount;
    self->elementCount = rowCount * columnCount;
    
    self->data = calloc(self->elementCount, sizeof(Float64));
    self->realData = self->data;
    self->imagData = calloc(self->elementCount, sizeof(Float64));
    self->tempData = calloc(self->elementCount, sizeof(Float64));
    self->tempImagData = calloc(self->elementCount, sizeof(Float64));
    
    self->label = "none";
    self->isComplex = &staticFalse;
    self->isSubmatrix = false;
    self->dataOnStack = false;
    memcpy((void *)&self->allocatedElementCount, &self->elementCount, sizeof(size_t));
    
    return self;
}

OVERLOADED Matrix *Matrix_new(size_t rowCount , size_t columnCount, char label[])
{
    Matrix *self = Matrix_new(rowCount, columnCount);
    self->label = label;
    
    return self;
}

OVERLOADED Matrix *Matrix_new(Matrix *input)
{
    
    Matrix *self = Matrix_new(input->rowCount, input->columnCount);
    
    Matrix_copy(input, self);
    
    return self;
}


OVERLOADED Matrix Matrix_new(Float64 data[], size_t rowCount, size_t columnCount)
{
    Matrix self = {0};
    
    self.elementCount = rowCount * columnCount;
    self.columnCount = columnCount;
    self.rowCount = rowCount;
    self.data = data;
    self.rowStride = columnCount;
    self.columnStride = rowCount;
    self.isComplex = &staticFalse;
    self.isSubmatrix = false;
    self.dataOnStack = true;
    memcpy((void *)&self.allocatedElementCount, &self.elementCount, sizeof(size_t));
    
    return self;
}

OVERLOADED Matrix Matrix_new(Float64 data[], size_t columnCount)
{
    Matrix self = Matrix_new(data, 1, columnCount);
    return self;
}

void Matrix_delete(Matrix *self)
{
    if (self->dataOnStack == false) {
        
        free(self->data);
        free(self->imagData);
        free(self->tempData);
        free(self->tempImagData);
        free(self);
    }
    
    *self = (Matrix){0};
}

void Matrix_scopedDelete(Matrix **self)
{
    Matrix_delete(*self);
}

#pragma mark - Setters / Getters -

OVERLOADED Float64 Matrix_getElement(Matrix *input, size_t index)
{
    _Matrix_checkForValidIndex(input, index);
    return input->data[index];
}

OVERLOADED Float64 Matrix_getElement(Matrix *input, RowOffset rowOffset, ColumnOffset columnOffset)
{
    return Matrix_getElement(input, input->rowStride * rowOffset.value + columnOffset.value);
}

OVERLOADED void Matrix_setElement(Matrix *input, size_t index, Float64 value)
{
    _Matrix_checkForValidIndex(input, index);
    input->data[index] = value;
}

OVERLOADED void Matrix_setElement(Matrix *input, RowOffset rowOffset, ColumnOffset columnOffset, Float64 value)
{
    Matrix_setElement(input, input->rowStride * rowOffset.value + columnOffset.value, value);
}

OVERLOADED void Matrix_setComplex(Matrix *input)
{
    input->isComplex = &staticTrue;
}

OVERLOADED void Matrix_setComplex(Matrix *input, bool state)
{
    if (state == true) {
        input->isComplex = &staticTrue;
    }
    else {
        input->isComplex = &staticFalse;
    }
}

#pragma mark - Safety -

OVERLOADED void Matrix_checkDimensionEquality(Matrix *inputA, Matrix *inputB)
{
    if (inputA->elementCount != inputB->elementCount
        ||
        inputA->rowCount != inputB->rowCount
        ||
        inputA->columnCount != inputB->columnCount) {
        
        printf("Matrix_checkDimensionEquality: Error, exiting\n");
        exit(-1);
    }
}

OVERLOADED void Matrix_checkDimensionEquality(size_t inputCount, ...)
{
    va_list argumentList;
    Matrix *inputArray[inputCount];
    va_start(argumentList, inputCount);
    
    for (size_t i = 0; i < inputCount; ++i) {
        
        inputArray[i] = va_arg(argumentList, Matrix *);
    }
    va_end(argumentList);
    
    for (size_t i = 0; i < inputCount - 1; ++i) {
        
        for (size_t j = i + 1; j < inputCount; ++j) {
            
            Matrix_checkDimensionEquality(inputArray[i], inputArray[j]);
        }
    }
}


void Matrix_checkElementCount(Matrix *self, size_t elementCount)
{
    if (self->elementCount < elementCount) {
        
        printf("Matrix_checkElementCount: Error, exiting\n");
        exit(-1);
    }
}

void Matrix_checkForVector(Matrix *self)
{
    if (self->rowCount != 1) {
        
        printf("Matrix_checkForVector: Error, exiting\n");
        exit(-1);
    }
}

void Matrix_checkSubmatrixDimensions(Matrix *input,
                                     size_t rowCount, size_t rowOffset,
                                     size_t columnCount, size_t columnOffset)
{
    if (rowCount + rowOffset > input->rowCount
        ||
        columnCount + columnOffset > input->columnCount) {
        
        printf("Matrix_checkSubmatrixDimensions: Error, exiting\n");
        exit(-1);
    }
}

Boolean Matrix_checkIsSubmatrix(Matrix *input)
{
    if (input->isSubmatrix == true) {
        
        return true;
    }
    else {
        
        return false;
    }
}

void Matrix_checkForValidIndex(Matrix *input, size_t index)
{
    if (index >= input->elementCount) {
        
        printf("Matrix_checkForValidIndex: Error, exiting\n");
        exit(-1);
    }
}

void Matrix_checkIsEmpty(Matrix *self)
{
    if (self->elementCount == 0) {
        
        printf("Matrix_checkIsEmpty: Matrix is empty, exiting\n");
        exit(-1);
    }
}

void Matrix_checkIsComplex(Matrix *self)
{
    if (*self->isComplex != true) {
        
        printf("Matrix_checkIsComplex: Error matrix is not complex, exiting");
        exit(-1);
    }
}

void Matrix_checkIsNotComplex(Matrix *self)
{
    if (*self->isComplex == true) {
        
        printf("Matrix_checkIsNotComplex: Error matrix is complex, exiting");
        exit(-1);
    }
}
void Matrix_checkCharacterEquality(char a, char b)
{
    if (a != b) {
        
        printf("Matrix_checkCharacterEquality: %c != %c, exiting\n", a, b);
        exit(-1);
    }
}

#pragma mark - Data copy -

OVERLOADED void Matrix_copy(Matrix source, Matrix destination)
{
    Matrix_copy(&source, &destination);
}

OVERLOADED void Matrix_copy(Matrix *source, Matrix *destination)
{
    _Matrix_checkDimensionEquality(2, source, destination);
    
    if (source->isSubmatrix == true
        ||
        destination->isSubmatrix == true) {
        
        for (size_t i = 0; i < source->rowCount; ++i) {
            
            cblas_dcopy((int)source->columnCount, Matrix_getRow(source, i), 1, Matrix_getRow(destination, i), 1);
        }
        
        if (Matrix_isComplex(source) == true) {
            
            for (size_t i = 0; i < source->rowCount; ++i) {
                
                cblas_dcopy((int)source->columnCount, Matrix_getImagRow(source, i), 1, Matrix_getImagRow(destination, i), 1);
            }
        }
        
    }
    else {
        
        cblas_dcopy((int)source->elementCount, source->data, 1, destination->data, 1);
        
        if (Matrix_isComplex(source) == true) {
            
            cblas_dcopy((int)source->elementCount, source->imagData, 1, destination->imagData, 1);
        }
    }
}


OVERLOADED void Matrix_copyRawData(Matrix *self, size_t elementCount, Float64 *data)
{
    _Matrix_checkElementCount(self, elementCount);
    
    if (self->isSubmatrix == true) {
        
        if (self->elementCount % elementCount == 0) {
            
            for (size_t i = 0; i < self->rowCount; ++i) {
                
                cblas_dcopy((int)self->columnCount, &data[self->columnCount * i], 1, Matrix_getRow(self, i), 1);
            }
        }
        else {
            
            size_t rowCount = elementCount / self->columnCount;
            size_t remainder = elementCount % (self->columnCount * rowCount);
            for (size_t i = 0; i < rowCount; ++i) {
                
                cblas_dcopy((int)self->columnCount, &data[self->columnCount * i], 1, Matrix_getRow(self, i), 1);
            }
            
            cblas_dcopy((int)remainder, &data[self->columnCount * rowCount], 1, Matrix_getRow(self, rowCount), 1);
        }
    }
    else {
        
        cblas_dcopy((int)elementCount, data, 1, self->data, 1);
    }
    
}

OVERLOADED void Matrix_copyRawComplexData(Matrix *self,
                                          size_t elementCount,
                                          Float64 *realData,
                                          Float64 *imagData)
{
    _Matrix_checkElementCount(self, elementCount);
    
    Matrix_isComplex(self) = true;
    
    if (self->isSubmatrix == true) {
        
        if (self->elementCount % elementCount == 0) {
            
            for (size_t i = 0; i < self->rowCount; ++i) {
                
                cblas_dcopy((int)self->columnCount, &realData[self->columnCount * i], 1, Matrix_getRow(self, i), 1);
                cblas_dcopy((int)self->columnCount, &imagData[self->columnCount * i], 1, Matrix_getImagRow(self, i), 1);
            }
        }
        else {
            
            size_t rowCount = elementCount / self->columnCount;
            size_t remainder = elementCount % (self->columnCount * rowCount);
            
            for (size_t i = 0; i < rowCount; ++i) {
                
                cblas_dcopy((int)self->columnCount, &realData[self->columnCount * i], 1, Matrix_getRow(self, i), 1);
                cblas_dcopy((int)self->columnCount, &imagData[self->columnCount * i], 1, Matrix_getImagRow(self, i), 1);
            }
            
            cblas_dcopy((int)remainder, &realData[self->columnCount * rowCount], 1, Matrix_getRow(self, rowCount), 1);
            cblas_dcopy((int)remainder, &imagData[self->columnCount * rowCount], 1, Matrix_getImagRow(self, rowCount), 1);
        }
    }
    else {
        
        cblas_dcopy((int)elementCount, realData, 1, self->data, 1);
        cblas_dcopy((int)elementCount, imagData, 1, self->imagData, 1);
    }
}

OVERLOADED void Matrix_subMatrixCopy(Matrix *source,
                                     RowCount rowCount,
                                     ColumnCount columnCount,
                                     Matrix *destination)
{
    Matrix_subMatrixCopy(source,
                         (SourceRowOffset){0},
                         (SourceColumnOffset){0},
                         rowCount,
                         columnCount,
                         destination,
                         (DestinationRowOffset){0},
                         (DestinationColumnOffset){0});
}

OVERLOADED void Matrix_subMatrixCopy(Matrix *source,
                                     ColumnCount columnCount,
                                     Matrix *destination)
{
    Matrix_subMatrixCopy(source,
                         (SourceRowOffset){0},
                         (SourceColumnOffset){0},
                         (RowCount){source->rowCount},
                         columnCount,
                         destination,
                         (DestinationRowOffset){0},
                         (DestinationColumnOffset){0});
}

OVERLOADED void Matrix_subMatrixCopy(Matrix *source,
                                     SourceRowOffset sourceStartRow,
                                     SourceColumnOffset sourceStartColumn,
                                     RowCount rowCount,
                                     ColumnCount columnCount,
                                     Matrix *destination,
                                     DestinationRowOffset destinationStartRow,
                                     DestinationColumnOffset destinationStartColumn)
{
    
    if (source->isSubmatrix == true || destination->isSubmatrix == true) {
        
        printf("Matrix_submatrixCopy: Error, source or destination cannot be a submatrix, exiting");
        exit(-1);
    }
    
    if (sourceStartRow.value >= source->rowCount
        ||
        sourceStartColumn.value >= source->columnCount
        ||
        sourceStartRow.value + rowCount.value > source->rowCount
        ||
        sourceStartColumn.value + columnCount.value > source->columnCount
        ||
        destinationStartRow.value >= destination->rowCount
        ||
        destinationStartColumn.value >= destination->columnCount
        ||
        destinationStartRow.value + rowCount.value > destination->rowCount
        ||
        destinationStartColumn.value + columnCount.value > destination->columnCount) {
        
        printf("Matrix_subMatrixCopy: Specified range illegal\nExiting\n");
        exit(-1);
    }
    
    vDSP_mmovD(&source->data[sourceStartColumn.value + (source->columnCount * sourceStartRow.value)],
               &destination->data[destinationStartColumn.value + (destination->columnCount * destinationStartRow.value)],
               columnCount.value,
               rowCount.value,
               source->columnCount,
               destination->columnCount);
    
    if (*source->isComplex == true) {
        
        vDSP_mmovD(&source->imagData[sourceStartColumn.value + (source->columnCount * sourceStartRow.value)],
                   &destination->imagData[destinationStartColumn.value + (destination->columnCount * destinationStartRow.value)],
                   columnCount.value,
                   rowCount.value,
                   source->columnCount,
                   destination->columnCount);
        Matrix_setComplex(destination);
    }
}

#pragma mark - Submatrix View -

Matrix Matrix_rowView(Matrix *input,
                      size_t rowOffset)
{
    return Matrix_submatrixView(input,
                                (RowCount){1},
                                (RowOffset){rowOffset},
                                (ColumnCount){input->columnCount},
                                (ColumnOffset){0});
}

Matrix Matrix_columnView(Matrix *input,
                         size_t columnOffset)
{
    return Matrix_submatrixView(input,
                                (RowCount){input->rowCount},
                                (RowOffset){0},
                                (ColumnCount){1},
                                (ColumnOffset){columnOffset});
}


OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                       RowCount rowCount)
{
    return Matrix_submatrixView(input,
                                rowCount,
                                (RowOffset){0},
                                (ColumnCount){input->columnCount},
                                (ColumnOffset){0});
}

OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                       RowCount rowCount,
                                       RowOffset rowOffset)
{
    return Matrix_submatrixView(input,
                                rowCount,
                                rowOffset,
                                (ColumnCount){input->columnCount},
                                (ColumnOffset){0});
}

OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                       ColumnCount columnCount)
{
    Matrix returnMatrix = Matrix_submatrixView(input,
                                               (RowCount){input->rowCount},
                                               (RowOffset){0},
                                               columnCount,
                                               (ColumnOffset){0});
    
    return returnMatrix;
}

OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                       RowCount rowCount,
                                       ColumnCount columnCount)
{
    return Matrix_submatrixView(input,
                                rowCount,
                                (RowOffset){0},
                                columnCount,
                                (ColumnOffset){0});
}

OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                       ColumnCount columnCount,
                                       ColumnOffset columnOffset)
{
    Matrix returnMatrix =  Matrix_submatrixView(input,
                                                (RowCount){input->rowCount},
                                                (RowOffset){0},
                                                columnCount,
                                                columnOffset);
    return returnMatrix;
}

OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                       ColumnOffset columnOffset)
{
    return Matrix_submatrixView(input,
                                (RowCount){input->rowCount},
                                (RowOffset){0},
                                (ColumnCount){input->columnCount - columnOffset.value},
                                columnOffset);
}

OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                       RowCount rowCount,
                                       ColumnCount columnCount,
                                       ColumnOffset columnOffset)
{
    return Matrix_submatrixView(input,
                                rowCount,
                                (RowOffset){0},
                                columnCount,
                                columnOffset);
}

OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                       RowCount rowCount,
                                       RowOffset rowOffset,
                                       ColumnCount columnCount)
{
    return Matrix_submatrixView(input,
                                rowCount,
                                rowOffset,
                                columnCount,
                                (ColumnOffset){0});
}

OVERLOADED Matrix Matrix_submatrixView(Matrix *input,
                                       RowCount rowCount, RowOffset rowOffset,
                                       ColumnCount columnCount, ColumnOffset columnOffset)
{
    _Matrix_checkSubmatrixDimensions(input, rowCount.value, rowOffset.value, columnCount.value, columnOffset.value);
    
    Matrix output = {0};
    
    output.rowCount = rowCount.value;
    output.columnCount = columnCount.value;
    output.elementCount = rowCount.value * columnCount.value;
    output.data = &input->data[rowOffset.value * input->rowStride + columnOffset.value];
    output.realData = output.data;
    output.imagData = &input->imagData[rowOffset.value * input->rowStride + columnOffset.value];
    output.tempData = &input->tempData[rowOffset.value * input->rowStride + columnOffset.value];
    output.tempImagData = &input->tempImagData[rowOffset.value * input->rowStride + columnOffset.value];
    output.rowStride = input->rowStride;
    output.columnStride = input->columnStride;
    output.isComplex = input->isComplex;
    output.isSubmatrix = true;
    size_t allocatedElementCount = 0;
    memcpy((void *)&output.allocatedElementCount, &allocatedElementCount, sizeof(size_t));
    
    return output;
}

OVERLOADED Matrix Matrix_submatrixView(Matrix *realPartMatrix, bool useRealPartForReal,
                                       Matrix *imagPartMatrix, bool useRealPartForImag,
                                       size_t rowCount, size_t rowOffset,
                                       size_t columnCount, size_t columnOffset)
{
    
    if (realPartMatrix->rowStride != imagPartMatrix->rowStride
        ||
        realPartMatrix->columnStride != imagPartMatrix->columnStride) {
        
        printf("Matrix_submatrixView: matrix strides are not equal, exiting");
        exit(-1);
    }
    
    _Matrix_checkDimensionEquality(2, realPartMatrix, imagPartMatrix);
    _Matrix_checkSubmatrixDimensions(realPartMatrix, rowCount, rowOffset, columnCount, columnOffset);
    
    Float64 *realPart;
    Float64 *imagPart;
    
    if (useRealPartForReal == false) {
        
        _Matrix_checkIsComplex(realPartMatrix);
        realPart = realPartMatrix->imagData;
    }
    else {
        
        realPart = realPartMatrix->data;
    }
    
    if (useRealPartForImag == false) {
        
        _Matrix_checkIsComplex(imagPartMatrix);
        imagPart = imagPartMatrix->imagData;
    }
    else {
        
        imagPart = imagPartMatrix->data;
    }
    Matrix output = {0};
    
    output.rowCount = rowCount;
    output.columnCount = columnCount;
    output.elementCount = rowCount * columnCount;
    output.data = &realPart[rowOffset * realPartMatrix->rowStride + columnOffset];
    output.imagData = &imagPart[rowOffset * imagPartMatrix->rowStride + columnOffset];
    output.rowStride = realPartMatrix->rowStride;
    output.columnStride = realPartMatrix->columnStride;
    output.isComplex = &staticTrue;
    output.isSubmatrix = true;
    size_t allocatedElementCount = 0;
    memcpy((void *)&output.allocatedElementCount, &allocatedElementCount, sizeof(size_t));
    
    return output;
}

Matrix Matrix_tempView(Matrix *input)
{
    Matrix output = {0};
    
    output.rowCount = input->rowCount;
    output.columnCount = input->columnCount;
    output.elementCount = input->elementCount;
    output.data = input->tempData;
    output.imagData = input->tempImagData;
    output.rowStride = input->rowStride;
    output.columnStride = input->columnStride;
    output.isComplex = input->isComplex;
    output.isSubmatrix = true;
    size_t allocatedElementCount = 0;
    memcpy((void *)&output.allocatedElementCount, &allocatedElementCount, sizeof(size_t));
    
    return output;
}

Matrix Matrix_realView(Matrix *input)
{
    Matrix output = {0};
    
    output.rowCount = input->rowCount;
    output.columnCount = input->columnCount;
    output.elementCount = input->elementCount;
    output.data = input->data;
    output.imagData = input->imagData;
    output.rowStride = input->rowStride;
    output.columnStride = input->columnStride;
    output.tempData = input->tempData;
    output.tempImagData = input->tempImagData;
    output.isComplex = &staticFalse;
    output.isSubmatrix = true;
    size_t allocatedElementCount = 0;
    memcpy((void *)&output.allocatedElementCount, &allocatedElementCount, sizeof(size_t));
    
    return output;
}

Matrix Matrix_imagView(Matrix *input)
{
    Matrix output = {0};
    
    output.rowCount = input->rowCount;
    output.columnCount = input->columnCount;
    output.elementCount = input->elementCount;
    output.data = input->imagData;
    output.imagData = input->data;
    output.rowStride = input->rowStride;
    output.columnStride = input->columnStride;
    output.tempData = input->tempData;
    output.tempImagData = input->tempImagData;
    output.isComplex = &staticFalse;
    output.isSubmatrix = true;
    size_t allocatedElementCount = 0;
    memcpy((void *)&output.allocatedElementCount, &allocatedElementCount, sizeof(size_t));
    
    return output;
}

Matrix Matrix_complexFlipView(Matrix *input)
{
    Matrix output = {0};
    
    output.rowCount = input->rowCount;
    output.columnCount = input->columnCount;
    output.elementCount = input->elementCount;
    output.data = input->imagData;
    output.imagData = input->data;
    output.rowStride = input->rowStride;
    output.columnStride = input->columnStride;
    output.tempData = input->tempData;
    output.tempImagData = input->tempImagData;
    output.isComplex = &staticTrue;
    output.isSubmatrix = true;
    size_t allocatedElementCount = 0;
    memcpy((void *)&output.allocatedElementCount, &allocatedElementCount, sizeof(size_t));
    
    return output;
}

void Matrix_complexFlip(Matrix **input)
{
    Float64 *temp;
    temp = (*input)->realData;
    (*input)->data = (*input)->realData = (*input)->imagData;
    (*input)->imagData = temp;
    temp = (*input)->tempData;
    (*input)->tempData = (*input)->tempImagData;
    (*input)->tempImagData = temp;
}
#pragma mark - Serialisation -

OVERLOADED void Matrix_saveHDF5File(Matrix self, const char *path, const char *dataSet)
{
    Matrix_saveHDF5File(&self, path, dataSet);
}

OVERLOADED void Matrix_saveHDF5File(Matrix *self, const char *path, const char *dataSet)
{
    if (self->isSubmatrix == true) {
        
        self = Matrix_new(self);
    }
    
    if(strncmp(path, "desktop", 7) == 0
       ||
       strncmp(path, "Desktop", 7) == 0) {
        
        char newPath[200];
        sprintf(newPath, "/Users/eddyc/Desktop%s.hdf", dataSet);
        path = newPath;
    }
    
    hid_t fileID, dataSetID, dataSpaceID;
    hsize_t dimensions[2];
    herr_t status;
    
    fileID = H5Fcreate(path,
                       H5F_ACC_TRUNC,
                       H5P_DEFAULT,
                       H5P_DEFAULT);
    
    dimensions[0] = self->rowCount;
    dimensions[1] = self->columnCount;
    
    dataSpaceID = H5Screate_simple(2, dimensions, NULL);
    
    if (*self->isComplex == false) {
        
        dataSetID = H5Dcreate(fileID, dataSet, H5T_NATIVE_DOUBLE, dataSpaceID, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        H5Dwrite (dataSetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, self->data);
        status = H5Dclose(dataSetID);
    }
    else {
        
        char dataSetReal[100]; char dataSetImag[100];
        sprintf(dataSetReal, "%s_real", dataSet);
        sprintf(dataSetImag, "%s_imag", dataSet);
        
        dataSetID = H5Dcreate(fileID, dataSetReal, H5T_NATIVE_DOUBLE, dataSpaceID, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite (dataSetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, self->data);
        status = H5Dclose(dataSetID);
        
        dataSetID = H5Dcreate(fileID, dataSetImag, H5T_NATIVE_DOUBLE, dataSpaceID, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite (dataSetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, self->imagData);
        status = H5Dclose(dataSetID);
    }
    
    status = H5Sclose(dataSpaceID);
    status = H5Fclose(fileID);
    
    if (self->isSubmatrix == true) {
        
        Matrix_delete(self);
    }
}


OVERLOADED Matrix *Matrix_openHDF5File(const char *path, const char *dataSet)
{
    if(strncmp(path, "desktop", 7) == 0
       ||
       strncmp(path, "Desktop", 7) == 0) {
        
        char newPath[200];
        sprintf(newPath, "/Users/eddyc/Desktop%s.hdf", dataSet);
        path = newPath;
    }
    
    hid_t fileID, dataSetID, dataSpaceID;
    hsize_t dimensions[2] = {0};
    
    herr_t status;
    
    fileID = H5Fopen(path,
                     H5F_ACC_RDONLY,
                     H5P_DEFAULT);
    
    dataSetID = H5Dopen2(fileID, dataSet, H5P_DEFAULT);
    dataSpaceID =  H5Dget_space(dataSetID);
    
    H5Sget_simple_extent_dims(dataSpaceID, dimensions, NULL);
    
    dimensions[0] = dimensions[0] == 0 ? 1 : dimensions[0];
    dimensions[1] = dimensions[1] == 0 ? 1 : dimensions[1];
    
    Matrix *self = Matrix_new((size_t)dimensions[0], (size_t)dimensions[1]);
    
    status = H5Dread(dataSetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, self->data);
    
    status = H5Sclose(dataSpaceID);
    status = H5Dclose(dataSetID);
    status = H5Fclose(fileID);
    
    return self;
}

void Matrix_addHDF5Dataset(Matrix *self, const char *path, const char *dataSet)
{
    
    hid_t fileID, dataSetID, dataSpaceID;
    hsize_t dimensions[2];
    herr_t status;
    
    fileID = H5Fopen(path,
                     H5F_ACC_RDWR,
                     H5P_DEFAULT);
    
    dimensions[0] = self->rowCount;
    dimensions[1] = self->columnCount;
    
    dataSpaceID = H5Screate_simple(2, dimensions, NULL);
    
    
    dataSetID = H5Dcreate(fileID,
                          dataSet,
                          H5T_NATIVE_DOUBLE,
                          dataSpaceID,
                          H5P_DEFAULT,
                          H5P_DEFAULT,
                          H5P_DEFAULT);
    
    H5Dwrite (dataSetID,
              H5T_NATIVE_DOUBLE,
              H5S_ALL,
              H5S_ALL,
              H5P_DEFAULT,
              self->data);
    
    status = H5Dclose(dataSetID);
    status = H5Sclose(dataSpaceID);
    status = H5Fclose(fileID);
    
}

#pragma mark - Printing -

OVERLOADED void Matrix_print(Matrix *input, CompareIndex comparison, size_t triggerIndex)
{
    if (comparison(triggerIndex) == true) {
        
        Matrix_print(input);
    }
}

OVERLOADED void Matrix_print(Matrix *input, const char *label)
{
    printf("\n%s\n", label);
    
    Matrix_print(input);
}

OVERLOADED void Matrix_print(Matrix input)
{
    Matrix_print(&input);
}

OVERLOADED void Matrix_print(Matrix *input)
{
    
    printf("---------------\n");
    printf("rows:%zd, columns:%zd\n", input->rowCount, input->columnCount);
    
    if (Matrix_isComplex(input) == false) {
        
        for (size_t i = 0; i < input->rowCount; ++i) {
          
            printf("[");
            for (size_t j = 0; j < input->columnCount; ++j) {
               
                if (j != input->columnCount - 1) {
                    
                    printf("%g, ", Matrix_getRow(input, i)[j]);
                }
                else {
                    
                    printf("%g", Matrix_getRow(input, i)[j]);
                }
            }
            
            printf("]\n");
        }
    }
    else {
        
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            for (size_t j = 0; j < input->columnCount; ++j) {
                
                Float64 real = Matrix_getRow(input, i)[j];
                Float64 imag = Matrix_getImagRow(input, i)[j];
                printf("[%g %c %gi]\t", real,
                       imag < 0 ? '-' : '+',
                       imag < 0 ? -imag : imag);
            }
            
            printf("\n");
        }
    }
    
    
    printf("---------------\n");
}

void Matrix_printRange(const Matrix *const self,
                       size_t startRow,
                       size_t rowCount,
                       size_t startColumn,
                       size_t columnCount)
{
    if (startRow + rowCount > self->rowCount
        ||
        startColumn + columnCount > self->columnCount) {
        
        printf("Matrix_printRange: Specified range illegal for matrix\nExiting\n");
        exit(-1);
    }
    
    printf("Matrix_printRange:\n\nrows = %zd, columns = %zd\n\nStart Row = %zd, Row Count = %zd, Start Column = %zd, Column Count = %zd\n", self->rowCount, self->columnCount, startRow, rowCount, startColumn, columnCount);
    
    if (Matrix_isComplex(self) == false) {
        
        for (size_t i = startRow; i < rowCount; ++i) {
            
            for (size_t j = startColumn; j < startColumn + columnCount; ++j) {
                
                printf("[%g]\t", self->data[i * self->columnCount + j]);
            }
            
            printf("\n");
        }
    }
    else {
        
        for (size_t i = startRow; i < rowCount; ++i) {
            
            for (size_t j = startColumn; j < startColumn + columnCount; ++j) {
                
                printf("[%g + %gi]", Matrix_getRow(self, i)[j], Matrix_getImagRow(self, i)[j]);
            }
            
            printf("\n");
        }
    }
    
    printf("\n\nEnd\n\n");
}

#pragma mark - Windowing functions -

OVERLOADED void Matrix_hanningWindow(Matrix *input)
{
    _Matrix_checkForVector(input);
    
    vDSP_hann_windowD(input->data, input->columnCount, vDSP_HANN_DENORM);
    
}

OVERLOADED void Matrix_sineWindow(Matrix *input)
{
    _Matrix_checkForVector(input);
    
    for (size_t i = 0; i < input->elementCount; ++i) {
        
        input->data[i] = sin((((Float64)i + 1.) - 0.5) * M_PI / input->elementCount);
    }
}

OVERLOADED void Matrix_sineWindow(Matrix *input, Float64 startValue, Float64 increment)
{
    _Matrix_checkForVector(input);
    
    for (size_t i = 0; i < input->columnCount; ++i) {
        
        input->data[i] = sin(startValue + i * increment);
    }
}

#pragma mark - Generate -

OVERLOADED void Matrix_reshape(Matrix *input, ColumnCount columns)
{
    Matrix_reshape(input, (RowCount){input->rowCount}, columns);
}

OVERLOADED void Matrix_reshape(Matrix *input, RowCount rows)
{
    Matrix_reshape(input, rows, (ColumnCount){input->columnCount});
}

OVERLOADED void Matrix_reshape(Matrix *input, RowCount rows, ColumnCount columns)
{
    size_t elementCount = rows.value * columns.value;
    
    if ((elementCount > input->allocatedElementCount
         &&
         input->isSubmatrix == false)
        ||
        (elementCount > input->elementCount
         &&
         input->isSubmatrix == true)) {
            
            printf("Matrix_reshape: Error, not enough allocated memory for specified size\n");
            exit(-1);
        }
    
    input->columnCount = columns.value;
    input->columnStride = rows.value;
    input->rowCount = rows.value;
    input->rowStride = columns.value;
    input->elementCount = elementCount;
}

OVERLOADED void Matrix_reshape(Matrix *input, Matrix *shape)
{
    Matrix_reshape(input, (RowCount){shape->rowCount}, (ColumnCount){shape->columnCount});
}

OVERLOADED void Matrix_transpose(Matrix *input)
{
    Matrix_transpose(input, input);
}

OVERLOADED void Matrix_transpose(Matrix *input, Matrix *output)
{
    if (input->isSubmatrix == true
        ||
        output->isSubmatrix == true) {
        
        printf("Matrix_transpose: Does not work with submatrices\n");
        exit(-1);
    }
    
    if (input->elementCount != output->elementCount) {
        
        printf("Matrix_transpose: element count not equal\n");
        exit(-1);
    }
    
    vDSP_mtransD(input->data, 1, input->tempData, 1, input->columnCount, input->rowCount);
    cblas_dcopy((UInt32)input->elementCount, input->tempData, 1, output->data, 1);
    
    size_t rowCount = input->rowCount;
    size_t columnCount = input->columnCount;
    size_t rowStride = input->rowStride;
    size_t columnStride = input->columnStride;
    output->rowCount= columnCount;
    output->columnCount = rowCount;
    output->rowStride = columnStride;
    output->columnStride = rowStride;
}

OVERLOADED void Matrix_fill(Matrix input, Float64 scalar)
{
    Matrix_fill(&input, scalar);
}

OVERLOADED void Matrix_fill(Matrix *input, Float64 scalar)
{
    if (input->isSubmatrix == true) {
        
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            vDSP_vfillD(&scalar, Matrix_getRow(input, i), 1, input->columnCount);
        }
    }
    else {
        
        vDSP_vfillD(&scalar, input->data, 1, input->elementCount);
    }
}

void Matrix_clear(const Matrix *const input)
{
    if (*input->isComplex == false) {
        
        Matrix_clearReal(input);
    }
    else {
        
        printf("Matrix_clear: error, complex clear not implemented yet");
        exit(-1);
    }
}
void Matrix_clearReal(const Matrix *const input)
{
    if (input->isSubmatrix == true) {
        
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            vDSP_vclrD(Matrix_getRow(input, i), 1, input->columnCount);
        }
    }
    else {
        
        vDSP_vclrD(input->data, 1, input->elementCount);
    }
}

OVERLOADED void Matrix_ramp(Matrix input, Float64 startValue, Float64 increment)
{
    Matrix_ramp(&input, startValue, increment);
}

OVERLOADED void Matrix_ramp(Matrix *input, Float64 startValue, Float64 increment)
{
    _Matrix_checkForVector(input);
    vDSP_vrampD(&startValue, &increment, input->data, 1, input->elementCount);
}


OVERLOADED Matrix *Matrix_ramp(Float64 startValue, Float64 endValue, Float64 increment)
{
    endValue += increment;
    size_t columnCount = (endValue - startValue) / increment;
    Matrix *result = Matrix_new(1, columnCount);
    vDSP_vrampD(&startValue, &increment, result->data, 1, result->elementCount);
    return result;
}

void Matrix_fillUsingIndexedElements(Matrix *input, Matrix *indices, Matrix *output)
{
    _Matrix_checkForVector(input);
    _Matrix_checkForVector(indices);
    
    for (size_t i = 0; i < indices->columnCount; ++i) {
        
        size_t index = indices->data[i];
        
        _Matrix_checkForValidIndex(input, index);
        
        output->data[i] = input->data[index];
    }
}

OVERLOADED void Matrix_fillElements(Matrix *input, CompareIndex comparison, Matrix *result)
{
    _Matrix_checkForVector(input);
    _Matrix_checkForVector(result);
    _Matrix_checkDimensionEquality(2, input, result);
    
    size_t resultIndex = 0;
    
    for (size_t i = 0; i < input->columnCount; ++i) {
        
        if (comparison(i) == true) {
            
            _Matrix_checkForValidIndex(result, resultIndex);
            result->data[resultIndex] = input->data[i];
            resultIndex++;
        }
    }
    
    result->columnCount = resultIndex;
    result->elementCount = resultIndex;
}


OVERLOADED void Matrix_fillElements(Matrix *input, CompareElement comparison, Matrix *result)
{
    _Matrix_checkForVector(input);
    _Matrix_checkForVector(result);
    
    size_t resultIndex = 0;
    
    for (size_t i = 0; i < input->columnCount; ++i) {
        
        if (comparison(input->data[i]) == true) {
            
            _Matrix_checkForValidIndex(result, resultIndex);
            result->data[resultIndex] = input->data[i];
            resultIndex++;
        }
    }
    
    result->columnCount = resultIndex;
    result->elementCount = resultIndex;
}

OVERLOADED void Matrix_fillElements(Matrix *input, CompareElement comparison, Matrix *comparisonData, Matrix *result)
{
    _Matrix_checkForVector(input);
    _Matrix_checkForVector(comparisonData);
    _Matrix_checkForVector(result);
    
    _Matrix_checkDimensionEquality(3, input, comparisonData, result);
    
    size_t resultIndex = 0;
    
    for (size_t i = 0; i < input->columnCount; ++i) {
        
        if (comparison(comparisonData->data[i]) == true) {
            
            _Matrix_checkForValidIndex(result, resultIndex);
            result->data[resultIndex] = input->data[i];
            resultIndex++;
        }
    }
    
    result->columnCount = resultIndex;
    result->elementCount = resultIndex;
}

OVERLOADED void Matrix_setElements(Matrix *input,
                                   CompareElement comparison,
                                   Float64 value,
                                   Matrix *result)
{
    _Matrix_checkForVector(input);
    _Matrix_checkForVector(result);
    _Matrix_checkDimensionEquality(2, input, result);
    
    for (size_t i = 0; i < input->columnCount; ++i) {
        
        if (comparison(input->data[i]) == true) {
            
            result->data[i] = value;
        }
    }
}

OVERLOADED void Matrix_setElements(Matrix *input,
                                   CompareElement comparison,
                                   Float64 (^value)(Float64 element),
                                   Matrix *result)
{
    _Matrix_checkForVector(input);
    _Matrix_checkForVector(result);
    _Matrix_checkDimensionEquality(2, input, result);
    
    for (size_t i = 0; i < input->columnCount; ++i) {
        
        if (comparison(input->data[i]) == true) {
            
            result->data[i] = value(result->data[i]);
        }
    }
}

void Matrix_fillIndices(Matrix *input, CompareElement comparison, Matrix *result)
{
    _Matrix_checkForVector(input);
    _Matrix_checkForVector(result);
    
    size_t resultIndex = 0;
    
    for (size_t i = 0; i < input->columnCount; ++i) {
        
        if (comparison(input->data[i]) == true) {
            
            _Matrix_checkForValidIndex(result, resultIndex);
            result->data[resultIndex] = i;
            resultIndex++;
        }
    }
    
    result->columnCount = resultIndex;
    result->elementCount = resultIndex;
}

void Matrix_insert(Matrix *input, Float64 value, size_t index)
{
    _Matrix_checkForVector(input);
    if (input->elementCount + 1 > input->allocatedElementCount) {
        
        printf("Matrix_insertValue error: not enough allocated memory\n");
        exit(-1);
    }
    
    const int size = (int)(input->columnCount - index);
    cblas_dcopy(size, &input->data[index], 1, input->tempData, 1);
    cblas_dcopy(size, input->tempData, 1, &input->data[index + 1], 1);
    
    input->elementCount++;
    input->columnCount++;
    input->rowStride++;
    
    input->data[index] = value;
    
}

OVERLOADED void Matrix_unique(Matrix *input)
{
    Matrix_unique(input, input);
}

OVERLOADED void Matrix_unique(Matrix *input, Matrix *output)
{
    _Matrix_checkForVector(input);
    _Matrix_checkDimensionEquality(2, input, output);
    _Matrix_checkIsNotComplex(input);
    _Matrix_checkIsNotComplex(output);
    
    cblas_dcopy((UInt32)input->columnCount, input->data, 1, input->tempData, 1);
    vDSP_vsortD(input->tempData, input->columnCount, 1);
    
    Float64 previousValue = NAN;
    size_t elementCount = 0;
    for (size_t i = 0; i < input->columnCount; ++i) {
        
        Float64 currentValue = input->tempData[i];
        
        if (previousValue == currentValue) {
            
            continue;
        }
        else {
            
            output->data[elementCount] = currentValue;
            previousValue = currentValue;
            elementCount++;
        }
    }
    
    output->elementCount = elementCount;
    output->columnCount = elementCount;
}

OVERLOADED void Matrix_complexConjugate(Matrix *input, Matrix *output)
{
    _Matrix_checkIsComplex(input);
    _Matrix_checkDimensionEquality(2, input, output);
    
    Matrix_setComplex(output);
    if (input->isSubmatrix == true
        ||
        output->isSubmatrix == true) {
        
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            cblas_dcopy((int)input->columnCount,
                        Matrix_getRow(input, i), 1,
                        Matrix_getRow(output, i), 1);
            
            vDSP_vnegD(Matrix_getImagRow(input, i), 1,
                       Matrix_getImagRow(output, i), 1,
                       input->columnCount);
        }
    }
    else {
        
        cblas_dcopy((int)input->columnCount,
                    input->data, 1,
                    output->data, 1);
        vDSP_vnegD(input->imagData, 1,
                   output->imagData, 1,
                   input->elementCount);
    }
}

const static size_t NSUM = 25;
double gaussrand()
{
    double x = 0;
    int i;
    for(i = 0; i < NSUM; i++)
        x += (double)rand() / RAND_MAX;
    
    x -= NSUM / 2.0;
    x /= sqrt(NSUM / 12.0);
    
    return x;
}

void Matrix_gaussian(Matrix *input, Float64 mean, Float64 variance)
{
    _Matrix_checkForVector(input);
    
    for (size_t i = 0; i < input->columnCount; ++i) {
        
        input->data[i] = gaussrand();
    }
    
    Matrix_multiply(input, variance);
    Matrix_add(input, mean);
}

void Matrix_diagonal(Matrix *input, Matrix *diagonal)
{
    _Matrix_checkForVector(diagonal);
    
    if (diagonal->elementCount != input->columnCount
        ||
        diagonal->elementCount != input->rowCount) {
        
        printf("Matrix_diagonal:Error, input dimensions not equal to diagonal\n Exiting...");
        exit(-1);
    }
    
    for (size_t i = 0; i < diagonal->elementCount; ++i) {
        
        Matrix_getRow(input, i)[i] = diagonal->data[i];
    }
}

void Matrix_identity(Matrix *input)
{
    if (input->rowCount != input->columnCount) {
        
        printf("Matrix_identity:Error, input dimensions not equal to diagonal\n Exiting...");
        exit(-1);
    }
    
    Matrix_clear(input);
    
    for (size_t i = 0; i < input->rowCount; ++i) {
        
        Matrix_getRow(input, i)[i] = 1;
    }
}

#pragma mark - Arithmetic functions -

OVERLOADED void Matrix_add(Matrix input, Float64 scalar)
{
    Matrix_add(&input, scalar);
}

OVERLOADED void Matrix_add(Matrix *input, Float64 scalar)
{
    Matrix_add(input, scalar, input);
}

OVERLOADED void Matrix_add(Matrix *input, Float64 scalar, Matrix *output)
{
    _Matrix_checkDimensionEquality(2, input, output);
    
    if (input->isSubmatrix == true
        ||
        output->isSubmatrix == true) {
        
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            vDSP_vsaddD(Matrix_getRow(input, i), 1,
                        &scalar,
                        Matrix_getRow(output, i), 1,
                        input->columnCount);
        }
    }
    else {
        
        vDSP_vsaddD(input->data, 1,
                    &scalar,
                    output->data, 1,
                    input->elementCount);
    }
}

OVERLOADED void Matrix_add(Matrix *inputA,
                           Matrix *inputB,
                           Matrix *result)
{
    _Matrix_checkDimensionEquality(3, inputA, inputB, result);
    
    if (inputA->isSubmatrix == true
        ||
        inputB->isSubmatrix == true
        ||
        result->isSubmatrix == true) {
        
        for (size_t i = 0; i < inputA->rowCount; ++i) {
            
            vDSP_vaddD(Matrix_getRow(inputA, i), 1, Matrix_getRow(inputB, i), 1,
                       Matrix_getRow(result, i), 1, inputA->columnCount);
        }
    }
    else {
        
        vDSP_vaddD(inputA->data, 1, inputB->data, 1,
                   result->data, 1, inputA->elementCount);
    }
}

OVERLOADED void Matrix_subtract(Matrix *input,
                                Matrix *negated,
                                Matrix *result)
{
    _Matrix_checkDimensionEquality(3, input, negated, result);
    
    if (input->isSubmatrix == true
        ||
        negated->isSubmatrix == true
        ||
        result->isSubmatrix == true) {
        
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            vDSP_vsubD(Matrix_getRow(input, i), 1, Matrix_getRow(negated, i), 1,
                       Matrix_getRow(result, i), 1, input->columnCount);
        }
    }
    else {
        
        vDSP_vsubD(input->data, 1, negated->data, 1,
                   result->data, 1, input->elementCount);
    }
}

OVERLOADED void Matrix_negate(Matrix *input)
{
    Matrix_negate(input, input);
}

OVERLOADED void Matrix_negate(Matrix *input, Matrix *output)
{
    _Matrix_checkDimensionEquality(2, input, output);
    
    if (*input->isComplex == true) {
        
        if (input->isSubmatrix == true
            ||
            output->isSubmatrix == true) {
            
            for (size_t i = 0; i < input->rowCount; ++i) {
                
                vDSP_vnegD(Matrix_getRow(input, i), 1,
                           Matrix_getRow(output, i), 1,
                           input->columnCount);
                vDSP_vnegD(Matrix_getImagRow(input, i), 1,
                           Matrix_getImagRow(output, i), 1,
                           input->columnCount);
            }
        }
        else {
            
            vDSP_vnegD(input->data, 1,
                       output->data, 1,
                       input->elementCount);
            vDSP_vnegD(input->imagData, 1,
                       output->imagData, 1,
                       input->elementCount);
        }
        
        Matrix_setComplex(output);
    }
    else {
        
        if (input->isSubmatrix == true
            ||
            output->isSubmatrix == true) {
            
            for (size_t i = 0; i < input->rowCount; ++i) {
                
                vDSP_vnegD(Matrix_getRow(input, i), 1,
                           Matrix_getRow(output, i), 1,
                           input->columnCount);
            }
        }
        else {
            
            vDSP_vnegD(input->data, 1,
                       output->data, 1,
                       input->elementCount);
        }
    }
}



OVERLOADED void Matrix_multiply(Matrix input, Float64 scalar)
{
    Matrix_multiply(&input, scalar);
}

OVERLOADED void Matrix_multiply(Matrix *input, Float64 scalar)
{
    Matrix_multiply(input, scalar, input);
}

OVERLOADED void Matrix_multiply(Matrix *input,
                                Float64 scalar,
                                Matrix *output)
{
    _Matrix_checkDimensionEquality(2, input, output);
    
    if (*input->isComplex == true) {
        
        if (input->isSubmatrix == true
            ||
            output->isSubmatrix == true) {
            
            for (size_t i = 0; i < input->rowCount; ++i) {
                
                vDSP_vsmulD(Matrix_getRow(input, i), 1,
                            &scalar,
                            Matrix_getRow(output, i), 1,
                            input->columnCount);
                vDSP_vsmulD(Matrix_getImagRow(input, i), 1,
                            &scalar,
                            Matrix_getImagRow(output, i), 1,
                            input->columnCount);
            }
        }
        else {
            
            vDSP_vsmulD(input->data, 1,
                        &scalar,
                        output->data, 1,
                        input->elementCount);
            vDSP_vsmulD(input->imagData, 1,
                        &scalar,
                        output->imagData, 1,
                        input->elementCount);
        }
        
        Matrix_setComplex(output);
    }
    else {
        
        if (input->isSubmatrix == true
            ||
            output->isSubmatrix == true) {
            
            for (size_t i = 0; i < input->rowCount; ++i) {
                
                vDSP_vsmulD(Matrix_getRow(input, i), 1,
                            &scalar,
                            Matrix_getRow(output, i), 1,
                            input->columnCount);
            }
        }
        else {
            
            vDSP_vsmulD(input->data, 1,
                        &scalar,
                        output->data, 1,
                        input->elementCount);
        }
    }
}

OVERLOADED void Matrix_multiply(Matrix *input,
                                double complex scalar)
{
    Matrix_multiply(input, scalar, input);
}

OVERLOADED void Matrix_multiply(Matrix *input,
                                double complex scalar,
                                Matrix *output)
{
    _Matrix_checkDimensionEquality(2, input, output);
    double realScalar = creal(scalar);
    double imagScalar = cimag(scalar);
    
    if (*input->isComplex == true) {
        
        if (input->isSubmatrix == true
            ||
            output->isSubmatrix == true) {
            
            for (size_t i = 0; i < input->rowCount; ++i) {
                
                vDSP_vsmulD(Matrix_getRow(input, i), 1,
                            &realScalar,
                            Matrix_getRow(output, i), 1,
                            input->columnCount);
                vDSP_vsmulD(Matrix_getImagRow(input, i), 1,
                            &imagScalar,
                            Matrix_getImagRow(output, i), 1,
                            input->columnCount);
            }
        }
        else {
            
            vDSP_vsmulD(input->data, 1,
                        &realScalar,
                        output->data, 1,
                        input->elementCount);
            vDSP_vsmulD(input->imagData, 1,
                        &imagScalar,
                        output->imagData, 1,
                        input->elementCount);
        }
    }
    else {
        
        if (input->isSubmatrix == true
            ||
            output->isSubmatrix == true) {
            
            for (size_t i = 0; i < input->rowCount; ++i) {
                
                vDSP_vsmulD(Matrix_getRow(input, i), 1,
                            &imagScalar,
                            Matrix_getImagRow(output, i), 1,
                            input->columnCount);
                vDSP_vsmulD(Matrix_getRow(input, i), 1,
                            &realScalar,
                            Matrix_getRow(output, i), 1,
                            input->columnCount);
            }
        }
        else {
            
            vDSP_vsmulD(input->data, 1,
                        &imagScalar,
                        output->imagData, 1,
                        input->elementCount);
            vDSP_vsmulD(input->data, 1,
                        &realScalar,
                        output->data, 1,
                        input->elementCount);
            
        }
    }
    
    Matrix_setComplex(output);
}


OVERLOADED void Matrix_multiply(Matrix *inputA,
                                Matrix *inputB,
                                Matrix *result)
{
    _Matrix_checkDimensionEquality(3, inputA, inputB, result);
    
    if (*inputA->isComplex == true || *inputB->isComplex == true) {
        
        if (*inputA->isComplex == true && *inputB->isComplex == true) {
            
            if (inputA->isSubmatrix == true
                ||
                inputB->isSubmatrix == true
                ||
                result->isSubmatrix == true) {
                
                for (size_t i = 0; i < inputA->rowCount; ++i) {
                    
                    DSPDoubleSplitComplex _inputA = {
                        .realp = Matrix_getRow(inputA, i),
                        .imagp = Matrix_getImagRow(inputA, i)
                    };
                    
                    DSPDoubleSplitComplex _inputB = {
                        .realp = Matrix_getRow(inputB, i),
                        .imagp = Matrix_getImagRow(inputB, i)
                    };
                    
                    DSPDoubleSplitComplex _result = {
                        .realp = Matrix_getRow(result, i),
                        .imagp = Matrix_getImagRow(result, i)
                    };
                    
                    vDSP_zvmulD(&_inputA, 1, &_inputB, 1, &_result, 1, result->elementCount, 1);
                }
            }
            else {
                
                DSPDoubleSplitComplex _inputA = {.realp = inputA->data, .imagp = inputA->imagData};
                DSPDoubleSplitComplex _inputB = {.realp = inputB->data, .imagp = inputB->imagData};
                DSPDoubleSplitComplex _result = {.realp = result->data, .imagp = result->imagData};
                
                vDSP_zvmulD(&_inputA, 1, &_inputB, 1, &_result, 1, result->elementCount, 1);
            }
        }
        else {
            
            if (*inputA->isComplex == true) {
                
                Matrix *temp = inputA;
                inputA = inputB;
                inputB = temp;
            }
            
            if (inputA->isSubmatrix == true
                ||
                inputB->isSubmatrix == true
                ||
                result->isSubmatrix == true) {
                
                for (size_t i = 0; i < inputA->rowCount; ++i) {
                    
                    vDSP_vmulD(Matrix_getRow(inputA, i), 1,
                               Matrix_getImagRow(inputB, i), 1,
                               Matrix_getImagRow(result, i), 1,
                               inputA->columnCount);
                    vDSP_vmulD(Matrix_getRow(inputA, i), 1,
                               Matrix_getRow(inputB, i), 1,
                               Matrix_getRow(result, i), 1,
                               inputA->columnCount);
                }
            }
            else {
                
                vDSP_vmulD(inputA->data, 1, inputB->imagData, 1, result->imagData, 1, result->elementCount);
                vDSP_vmulD(inputA->data, 1, inputB->data, 1, result->data, 1, result->elementCount);
            }
        }
        
        Matrix_setComplex(result);
    }
    else {
        
        if (inputA->isSubmatrix == true
            ||
            inputB->isSubmatrix == true
            ||
            result->isSubmatrix == true) {
            
            for (size_t i = 0; i < inputA->rowCount; ++i) {
                
                vDSP_vmulD(Matrix_getRow(inputA, i), 1,
                           Matrix_getRow(inputB, i), 1,
                           Matrix_getRow(result, i), 1,
                           inputA->columnCount);
            }
        }
        else {
            
            vDSP_vmulD(inputA->data, 1, inputB->data, 1, result->data, 1, result->elementCount);
        }
    }
}

void Matrix_dotProduct(Matrix *inputA, Matrix *inputB, Matrix *output)
{
    if (inputA->columnCount != inputB->rowCount
        ||
        inputB->columnCount != output->columnCount
        ||
        inputA->rowCount != output->rowCount) {
        
        printf("Matrix_dotProduct: Error, Exiting\n");
        exit(-1);
    }
    
    if (inputA == output || inputB == output) {
        
        printf("Matrix_dotProduct: Error, Exiting\n");
        exit(-1);
    }
    
    vDSP_mmulD(inputA->data, 1, inputB->data, 1, output->data, 1, inputA->rowCount, inputB->columnCount, inputA->columnCount);
}

OVERLOADED void Matrix_divide(Matrix input, Float64 scalar)
{
    Matrix_divide(&input, scalar);
}

OVERLOADED void Matrix_divide(Matrix *input, Float64 scalar)
{
    Matrix_divide(input, scalar, input);
}

OVERLOADED void Matrix_divide(Matrix *input,
                              Float64 scalar,
                              Matrix *output)
{
    _Matrix_checkDimensionEquality(2, input, output);
    
    if (input->isSubmatrix == true
        ||
        output->isSubmatrix == true) {
        
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            vDSP_vsdivD(Matrix_getRow(input, i), 1,
                        &scalar,
                        Matrix_getRow(output, i), 1,
                        input->columnCount);
        }
    }
    else {
        
        vDSP_vsdivD(input->data, 1,
                    &scalar,
                    output->data, 1,
                    input->elementCount);
    }
}

OVERLOADED void Matrix_divide(Matrix *inputA,
                              Matrix *inputB,
                              Matrix *result)
{
    _Matrix_checkDimensionEquality(3, inputA, inputB, result);
    
    if (inputA->isSubmatrix == true
        ||
        inputB->isSubmatrix == true
        ||
        result->isSubmatrix == true) {
        
        for (size_t i = 0; i < inputA->rowCount; ++i) {
            
            vDSP_vdivD(Matrix_getRow(inputA, i), 1, Matrix_getRow(inputB, i), 1, Matrix_getRow(result, i), 1, inputA->columnCount);
        }
    }
    else {
        
        vDSP_vdivD(inputA->data, 1, inputB->data, 1, result->data, 1, result->elementCount);
    }
}

OVERLOADED Float64 Matrix_sum(Matrix *input, size_t startIndex, size_t size)
{
    _Matrix_checkForVector(input);
    _Matrix_checkForValidIndex(input, startIndex + size - 1);
    
    Float64 output = 0;
    vDSP_sveD(&input->data[startIndex], 1, &output, size);
    
    return output;
}

OVERLOADED Float64 Matrix_sum(Matrix *input)
{
    Float64 output = Matrix_sum(input, 0, input->columnCount);
    
    return output;
}

OVERLOADED Float64 Matrix_sum(Matrix *input, Matrix *indices)
{
    _Matrix_checkForVector(input);
    _Matrix_checkForVector(indices);
    
    Float64 sum = 0;
    
    for (size_t i = 0; i < indices->columnCount; ++i) {
        
        size_t index = indices->data[i];
        
        _Matrix_checkForValidIndex(input, index);
        
        sum += input->data[index];
    }
    
    return sum;
}

void Matrix_sumRows(Matrix *input, Matrix *output)
{
    _Matrix_checkForVector(output);
    
    if (input->columnCount != output->columnCount) {
        
        printf("Matrix_sumRows: output matrix row count not equal to input matrix row count\n");
        exit(-1);
    }
    
    for (size_t i = 0; i < input->rowCount; ++i) {
       
        vDSP_vaddD(&input->data[input->columnCount * i], 1, output->data, 1, output->data, 1, output->columnCount);
    }

}

OVERLOADED void Matrix_difference(Matrix *input, Matrix *output)
{
    _Matrix_checkForVector(input);
    _Matrix_checkForVector(output);
    
    if (output->columnCount != input->columnCount - 1) {
        
        printf("Matrix_difference: output matrix needs to be one element smaller than input matrix\n");
        exit(-1);
    }
    
    vDSP_vsubD(input->data, 1, &input->data[1], 1, output->data, 1, input->columnCount - 1);
}

Matrix *Matrix_newUnique(Matrix *input)
{
    _Matrix_checkForVector(input);
    Matrix *self = Matrix_new(input->rowCount, input->columnCount);
    
    cblas_dcopy((UInt32)input->columnCount, input->data, 1, input->tempData, 1);
    vDSP_vsortD(input->tempData, input->columnCount, 1);
    
    Float64 previousValue = NAN;
    size_t elementCount = 0;
    for (size_t i = 0; i < input->columnCount; ++i) {
        
        Float64 currentValue = input->tempData[i];
        
        if (previousValue == currentValue) {
            
            continue;
        }
        else {
            
            self->data[elementCount] = currentValue;
            previousValue = currentValue;
            elementCount++;
        }
    }
    
    self->elementCount = elementCount;
    self->columnCount = elementCount;
    
    return self;
}

OVERLOADED Matrix *Matrix_newByConcatenatingMatrices(Matrix *start, Matrix *end)
{
    _Matrix_checkForVector(start);
    _Matrix_checkForVector(end);
    
    Matrix *self = Matrix_new(1, start->columnCount + end->columnCount);
    
    cblas_dcopy((UInt32)start->columnCount, start->data, 1, self->data, 1);
    cblas_dcopy((UInt32)end->columnCount, end->data, 1, &self->data[start->columnCount], 1);
    
    return self;
}

Matrix *Matrix_newDifference(Matrix *inputA, Matrix *inputB)
{
    _Matrix_checkForVector(inputA);
    _Matrix_checkForVector(inputB);
    
    _Matrix *uniqueA = Matrix_newUnique(inputA);
    _Matrix *uniqueB = Matrix_newUnique(inputB);
    
    _Matrix *input = Matrix_newByConcatenatingMatrices(uniqueA, uniqueB);
    
    Matrix *output = Matrix_new(1, input->columnCount);
    Matrix_sort(input, input, true);
    
    size_t columnCount = 0;
    
    for (size_t i = 1; i < input->columnCount; ++i) {
        
        size_t elementA = input->data[i - 1];
        size_t elementB = input->data[i];
        
        if (elementA == elementB) {
            
            i++;
            continue;
        }
        else {
            
            output->data[columnCount] = input->data[i - 1];
            columnCount++;
        }
    }
    
    size_t elementA = input->data[input->columnCount - 2];
    size_t elementB = input->data[input->columnCount - 1];
    
    if (elementA != elementB) {
        
        output->data[columnCount] = input->data[input->columnCount - 1];
        columnCount++;
    }
    
    output->columnCount = columnCount;
    
    return output;
}

OVERLOADED void Matrix_difference(Matrix *inputA, Matrix *inputB, Matrix *output)
{
    _Matrix_checkForVector(inputA);
    _Matrix_checkForVector(inputB);
    
    Matrix_reshape(output, inputA);
    Matrix_unique(inputA, output);
    
    
    size_t aSize = output->columnCount;
    Matrix_reshape(output, (RowCount){1}, (ColumnCount){output->allocatedElementCount});
    Matrix outputView = Matrix_submatrixView(output, (RowCount){1}, (ColumnCount){inputB->columnCount}, (ColumnOffset){aSize});
    Matrix_unique(inputB, &outputView);
    
    Matrix_reshape(output, (RowCount){1}, (ColumnCount){aSize + outputView.columnCount});
    Matrix_sort(output, output, true);
    
    size_t columnCount = 0;
    
    for (size_t i = 1; i < output->columnCount; ++i) {
        
        Float64 elementA = output->data[i - 1];
        Float64 elementB = output->data[i];
        
        if (elementA == elementB) {
            
            i++;
            continue;
        }
        else {
            
            output->tempData[columnCount] = output->data[i - 1];
            columnCount++;
        }
    }
    
    Float64 elementA = output->data[output->columnCount - 2];
    Float64 elementB = output->data[output->columnCount - 1];
    
    if (elementA != elementB) {
        
        output->tempData[columnCount] = output->data[output->columnCount - 1];
        columnCount++;
    }
    
    output->columnCount = columnCount;
    
    cblas_dcopy((int)columnCount, output->tempData, 1, output->data, 1);
    
}

OVERLOADED void Matrix_floor(Matrix *input)
{
    Matrix_floor(input, input);
}

OVERLOADED void Matrix_floor(Matrix *input, Matrix *output)
{
    _Matrix_checkDimensionEquality(2, input, output);
    
    if (input->isSubmatrix == true
        ||
        output->isSubmatrix == true) {
        
        int columnCount = (int)input->columnCount;
        
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            vvfloor(Matrix_getRow(output, i), Matrix_getRow(input, i), &columnCount);
        }
    }
    else {
        
        int elementCount = (int)input->elementCount;
        vvfloor(output->data, input->data, &elementCount);
    }
}

OVERLOADED void Matrix_ceiling(Matrix *input)
{
    Matrix_ceiling(input, input);
}

OVERLOADED void Matrix_cummulativeSum(Matrix *input)
{
    Matrix_cummulativeSum(input, input);
}

OVERLOADED void Matrix_cummulativeSum(Matrix *input, Matrix *output)
{
    _Matrix_checkForVector(input);
    _Matrix_checkDimensionEquality(2, input, output);
    
    Float64 one = 1;
    Float64 firstElement = input->data[0];
    vDSP_vrsumD(input->data, 1, &one, output->data, 1, input->columnCount);
    vDSP_vsaddD(output->data, 1, &firstElement, output->data, 1, output->columnCount);
}

OVERLOADED void Matrix_ceiling(Matrix *input, Matrix *output)
{
    _Matrix_checkDimensionEquality(2, input, output);
    
    if (input->isSubmatrix == true
        ||
        output->isSubmatrix == true) {
        
        int columnCount = (int)input->columnCount;
        
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            vvceil(Matrix_getRow(output, i), Matrix_getRow(input, i), &columnCount);
        }
    }
    else {
        
        int elementCount = (int)input->elementCount;
        vvceil(output->data, input->data, &elementCount);
    }
}


#pragma mark - Algebraic functions -

OVERLOADED void Matrix_square(Matrix *input)
{
    Matrix_square(input, input);
}

OVERLOADED void Matrix_square(Matrix *input, Matrix *output)
{
    _Matrix_checkDimensionEquality(2, input, output);
    
    if (input->isSubmatrix == true
        ||
        output->isSubmatrix == true) {
        
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            vDSP_vsqD(Matrix_getRow(input, i), 1, Matrix_getRow(output, i), 1, input->columnCount);
        }
    }
    else {
        
        vDSP_vsqD(input->data, 1, output->data, 1, input->elementCount);
    }
}

Float64 Matrix_sumOfSquares(Matrix *input)
{
    Float64 totalSum = 0;
    
    if (input->isSubmatrix == true) {
        
        Float64 currentSum = 0;
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            vDSP_svesqD(Matrix_getRow(input, i), 1, &currentSum, input->columnCount);
            totalSum += currentSum;
        }
    }
    else {
        
        vDSP_svesqD(input->data, 1, &totalSum, input->elementCount);
    }
    return totalSum;
}

OVERLOADED void Matrix_squareRoot(Matrix *input)
{
    Matrix_squareRoot(input, input);
}

OVERLOADED void Matrix_squareRoot(Matrix *input, Matrix *output)
{
    _Matrix_checkDimensionEquality(2, input, output);
    
    if (Matrix_isComplex(input) == false) {
        
        if (input->isSubmatrix == true
            ||
            output->isSubmatrix == true) {
            
            for (size_t i = 0; i < input->rowCount; ++i) {
                
                int columnCount = (int)input->columnCount;
                vvsqrt(Matrix_getRow(output, i), Matrix_getRow(input, i), &columnCount);
            }
        }
        else {
            
            int elementCount = (int)input->elementCount;
            vvsqrt(output->data, input->data, &elementCount);
        }
    }
    else {
        //TODO: Complex square root
        printf("Complex square root not implemented\n");
        exit(-1);
    }
    
}


void Matrix_modulus(Matrix *input, Matrix *output)
{
    _Matrix_checkDimensionEquality(2, input, output);
    
    if (Matrix_isComplex(input) == false) {
        
        if (input->isSubmatrix == true
            ||
            output->isSubmatrix == true) {
            
            for (size_t i = 0; i < input->rowCount; ++i) {
                
                vDSP_vabsD(Matrix_getRow(input, i), 1, Matrix_getRow(output, i), 1, input->columnCount);
            }
        }
        else {
            
            vDSP_vabsD(input->data, 1, output->data, 1, input->elementCount);
        }
    }
    else {
        
        DSPDoubleSplitComplex inputComplex;
        
        if (input->isSubmatrix == true
            ||
            output->isSubmatrix == true) {
            
            for (size_t i = 0; i < input->rowCount; ++i) {
                
                inputComplex.realp = Matrix_getRow(input, i);
                inputComplex.imagp = Matrix_getImagRow(input, i);
                
                vDSP_vdistD(Matrix_getRow(input, i), 1, Matrix_getImagRow(input, i), 1, output->data, 1, input->columnCount);
            }
        }
        else {
            
            inputComplex.realp = input->data;
            inputComplex.imagp = input->imagData;
            vDSP_vdistD(input->data, 1, input->imagData, 1, output->data, 1, input->elementCount);
        }
    }
    
    Matrix_isComplex(output) = false;
}


OVERLOADED void Matrix_power(Matrix *input, Matrix *power, Matrix *output)
{
    _Matrix_checkDimensionEquality(3, input, power, output);
    
    if (input->isSubmatrix
        ||
        power->isSubmatrix
        ||
        output->isSubmatrix) {
        
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            int columnCount = (int)input->columnCount;
            vvpow(Matrix_getRow(output, i), Matrix_getRow(power, i), Matrix_getRow(input, i), &columnCount);
        }
    }
    else {
        
        int count = (int)input->elementCount;
        vvpow(output->data, power->data, input->data, &count);
    }
}

OVERLOADED void Matrix_log10(Matrix *input)
{
    Matrix_log10(input, input);
}

OVERLOADED void Matrix_log10(Matrix *input, Matrix *output)
{
    _Matrix_checkDimensionEquality(2, input, output);
    
    if (input->isSubmatrix
        ||
        output->isSubmatrix) {
        
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            int columnCount = (int)input->columnCount;
            vvlog10(Matrix_getRow(output, i), Matrix_getRow(input, i), &columnCount);
        }
    }
    else {
        
        int count = (int)input->elementCount;
        vvlog10(output->data, input->data, &count);
    }
}

void Matrix_exponent(Matrix *input, Matrix *output)
{
    _Matrix_checkDimensionEquality(2, input, output);
    if (input->isSubmatrix
        ||
        output->isSubmatrix) {
        
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            int columnCount = (int)input->columnCount;
            vvexp(Matrix_getRow(output, i), Matrix_getRow(input, i), &columnCount);
        }
    }
    else {
        
        int count = (int)input->elementCount;
        vvexp(output->data, input->data, &count);
    }
    
}
void Matrix_complexExponent(Matrix *input, size_t frameNumber)
{
    _Matrix_checkIsComplex(input);
    
    Matrix temp = Matrix_tempView(input);
    Matrix inputImag = Matrix_imagView(input);
    Matrix inputReal = Matrix_realView(input);
    Matrix tempImag = Matrix_imagView(&temp);
    
    Matrix_cosine(&inputImag, &temp);
    Matrix_sine(&inputImag, &tempImag);
    Matrix_exponent(input, input);
    
    Matrix_multiply(&inputReal, &temp, input);
}


void Matrix_complexExponentBlankReal(Matrix *input)
{
    
    Matrix inputImag = Matrix_imagView(input);
    Matrix inputReal = Matrix_realView(input);
    
    Matrix_cosine(&inputImag, &inputReal);
    Matrix_sine(&inputImag, &inputImag);
}
#pragma mark - Trigonometric functions -

void Matrix_angle(Matrix *input, Matrix *result)
{
    _Matrix_checkForVector(input);
    _Matrix_checkIsComplex(input);
    _Matrix_checkDimensionEquality(2, input, result);
    
    DSPDoubleSplitComplex inputComplex = {.realp = input->data, .imagp = input->imagData};
    vDSP_zvphasD(&inputComplex, 1, result->data, 1, result->elementCount);
}

void Matrix_sine(Matrix *input, Matrix *result)
{
    
    _Matrix_checkDimensionEquality(2, input, result);
    
    if (input->isSubmatrix
        ||
        result->isSubmatrix) {
        
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            int columnCount = (int)input->columnCount;
            vvsin(Matrix_getRow(result, i), Matrix_getRow(input, i), &columnCount);
        }
    }
    else {
        
        int count = (int)input->elementCount;
        vvsin(result->data, input->data, &count);
    }
}

void Matrix_cosine(Matrix *input, Matrix *result)
{
    
    _Matrix_checkDimensionEquality(2, input, result);
    
    if (input->isSubmatrix
        ||
        result->isSubmatrix) {
        
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            int columnCount = (int)input->columnCount;
            vvcos(Matrix_getRow(result, i), Matrix_getRow(input, i), &columnCount);
        }
    }
    else {
        
        int count = (int)input->elementCount;
        vvcos(result->data, input->data, &count);
    }
    
}

void Matrix_magnitudeSquared(Matrix *input, Matrix *output)
{
    _Matrix_checkForVector(input);
    _Matrix_checkIsComplex(input);
    _Matrix_checkDimensionEquality(2, input, output);
    
    DSPDoubleSplitComplex tempComplex = {.realp = input->data, .imagp = input->imagData};
    vDSP_zvmagsD(&tempComplex, 1, output->data, 1, output->columnCount);
}

#pragma mark - Information -

void Matrix_minimum(Matrix *input, Float64 *minimum, size_t *index)
{
    _Matrix_checkIsEmpty(input);
    _Matrix_checkForVector(input);
    
    vDSP_minviD(input->data, 1, minimum, index, input->elementCount);
}

void Matrix_maximum(Matrix *input, Float64 *maximum, size_t *index)
{
    _Matrix_checkForVector(input);
    _Matrix_checkIsEmpty(input);
    
    vDSP_maxviD(input->data, 1, maximum, index, input->elementCount);
}

void Matrix_domain(Matrix *input, Float64 domain[2])
{
    size_t index;
    Matrix_minimum(input, &domain[0], &index);
    Matrix_maximum(input, &domain[1], &index);
}

size_t Matrix_countElements(Matrix *input, bool (^comparison)(Float64))
{
    size_t elementCount = 0;
    for (size_t i = 0; i < input->elementCount; ++i) {
        
        if (comparison(input->data[i]) == true) {
            
            elementCount++;
        }
    }
    
    return elementCount;
}

#pragma mark - Other functions -

OVERLOADED void Matrix_interpolate(Matrix *input, Matrix *indices, Matrix *output)
{
    _Matrix_checkDimensionEquality(2, indices, output);
    _Matrix_checkForVector(indices);
    
    for (size_t i = 0; i < input->rowCount; ++i) {
        
        Float64 *currentInputRow = Matrix_getRow(input, i);
        Float64 *currentOutputRow = Matrix_getRow(output, i);
        
        vDSP_vlintD(currentInputRow, indices->data, 1, currentOutputRow, 1, output->columnCount, input->columnCount);
    }
}

Float64 getInterpolatedValue(Float64 startIndex, Float64 startValue, Float64 endIndex, Float64 endValue, Float64 inputIndex)
{
    if (startIndex < 0
        ||
        startIndex > endIndex
        ||
        inputIndex < startIndex
        ||
        inputIndex > endIndex) {
        
        printf("getInterpolatedValue: Illegal arguments\n");
        
        return NAN;
    }
    Float64 indexDifference = endIndex - startIndex;
    Float64 valueDifference = endValue - startValue;
    
    Float64 inputIndexPercent = (inputIndex - startIndex) / indexDifference;
    Float64 result = startValue + (valueDifference * inputIndexPercent);
    
    return result;
}

OVERLOADED void Matrix_interpolate(Matrix *currentIndices, Matrix *currentValues,
                                   Matrix *sortIndices, Matrix *sortedCurrentIndices, Matrix *sortedCurrentValues,
                                   Matrix *newIndices, Matrix *output)
{
    
    
    Float64 nan = 0;
    _Matrix_checkDimensionEquality(7, currentIndices, currentValues, sortIndices,
                                   sortedCurrentIndices, sortedCurrentValues, newIndices, output);
    
    Matrix_getSortIndices(currentIndices, sortIndices, true);
    Matrix_fillUsingIndexedElements(currentIndices, sortIndices, sortedCurrentIndices);
    Matrix_fillUsingIndexedElements(currentValues, sortIndices, sortedCurrentValues);
    
    Float64 newIndicesMinimum;
    Float64 newIndicesMaximum;
    
    Float64 currentIndicesMinimum;
    Float64 currentIndicesMaximum;
    
    size_t index = 0;
    Matrix_minimum(sortedCurrentIndices, &currentIndicesMinimum, &index);
    Matrix_maximum(sortedCurrentIndices, &currentIndicesMaximum, &index);
    Matrix_minimum(newIndices, &newIndicesMinimum, &index);
    Matrix_maximum(newIndices, &newIndicesMaximum, &index);
    
    if (newIndicesMinimum > currentIndicesMaximum
        ||
        newIndicesMaximum < currentIndicesMinimum) {
        
        Matrix_fill(output, nan);
        return;
    }
    
    for (size_t row = 0; row < output->rowCount; ++row) {
        
        Float64 *currentOutputRow = Matrix_getRow(output, row);
        
        for (size_t i = 0, valuePosition = 0; i < newIndices->columnCount; ++i) {
            
            Float64 currentNewIndex = newIndices->data[i];
            
            if (currentNewIndex < currentIndicesMinimum
                ||
                currentNewIndex > currentIndicesMaximum) {
                
                currentOutputRow[i] = nan;
            }
            else {
                
                for (; valuePosition < sortedCurrentIndices->columnCount - 1; ++valuePosition) {
                    
                    Float64 currentStartIndex; Float64 currentStartValue;
                    Float64 currentEndIndex; Float64 currentEndValue;
                    
                    currentStartIndex = sortedCurrentIndices->data[valuePosition];
                    currentStartValue = sortedCurrentValues->data[valuePosition];
                    
                    currentEndIndex = sortedCurrentIndices->data[valuePosition + 1];
                    currentEndValue = sortedCurrentValues->data[valuePosition + 1];
                    
                    if (currentNewIndex >= currentStartIndex
                        &&
                        currentNewIndex <= currentEndIndex) {
                        
                        currentOutputRow[i] = getInterpolatedValue(currentStartIndex, currentStartValue, currentEndIndex, currentEndValue, currentNewIndex);
                        break;
                    }
                }
            }
        }
    }
}

OVERLOADED void Matrix_fastInterpolate(Matrix *currentIndices, Matrix *currentValues,
                                       Matrix *newIndices, Matrix *output)
{
    
    Float64 nan = 0;
    
    Float64 newIndicesMinimum;
    Float64 newIndicesMaximum;
    
    Float64 currentIndicesMinimum;
    Float64 currentIndicesMaximum;
    
    size_t index = 0;
    Matrix_minimum(currentIndices, &currentIndicesMinimum, &index);
    Matrix_maximum(currentIndices, &currentIndicesMaximum, &index);
    Matrix_minimum(newIndices, &newIndicesMinimum, &index);
    Matrix_maximum(newIndices, &newIndicesMaximum, &index);
    
    if (newIndicesMinimum > currentIndicesMaximum
        ||
        newIndicesMaximum < currentIndicesMinimum) {
        
        Matrix_fill(output, nan);
        return;
    }
    
    for (size_t row = 0; row < output->rowCount; ++row) {
        
        Float64 *currentOutputRow = Matrix_getRow(output, row);
        Float64 *currentValuesRow = Matrix_getRow(currentValues, row);
        
        for (size_t i = 0, valuePosition = 0; i < newIndices->columnCount; ++i) {
            
            Float64 currentNewIndex = newIndices->data[i];
            
            if (currentNewIndex < currentIndicesMinimum
                ||
                currentNewIndex > currentIndicesMaximum) {
                
                currentOutputRow[i] = nan;
            }
            else {
                
                for (; valuePosition < currentIndices->columnCount - 2; ++valuePosition) {
                    
                    Float64 currentStartIndex; Float64 currentStartValue;
                    Float64 currentEndIndex; Float64 currentEndValue;
                    
                    currentStartIndex = currentIndices->data[valuePosition];
                    currentStartValue = currentValuesRow[valuePosition];
                    
                    currentEndIndex = currentIndices->data[valuePosition + 1];
                    currentEndValue = currentValuesRow[valuePosition + 1];
                    
                    if (currentNewIndex >= currentStartIndex
                        &&
                        currentNewIndex <= currentEndIndex) {
                        
                        currentOutputRow[i] = getInterpolatedValue(currentStartIndex, currentStartValue, currentEndIndex, currentEndValue, currentNewIndex);
                        break;
                    }
                }
            }
        }
    }
}

OVERLOADED void Matrix_getSortIndices(Matrix *input, Matrix *indices, bool forward)
{
    if (input->isSubmatrix
        ||
        indices->isSubmatrix) {
        
        printf("Matrix_getSortIndices: Does not work with submatrices, exiting\n");
        exit(-1);
    }
    
    _Matrix_checkDimensionEquality(2, input, indices);
    
    int order = forward == true ? 1 : -1;
    
    size_t temp[indices->elementCount];
    for (size_t i = 0; i < indices->elementCount; ++i) {
        
        temp[i] = i;
    }
    
    vDSP_vsortiD(input->data, temp, NULL, input->elementCount, order);
    
    for (size_t i = 0; i < indices->elementCount; ++i) {
        
        indices->data[i] = temp[i];
    }
}

OVERLOADED void Matrix_sort(Matrix *input, Boolean forward)
{
    Matrix_sort(input, input, forward);
}

OVERLOADED void Matrix_sort(Matrix *input, Matrix *output, Boolean forward)
{
    _Matrix_checkForVector(input);
    _Matrix_checkDimensionEquality(2, input, output);
    
    int order = forward == true ? 1 : -1;
    
    if (input == output) {
        
        vDSP_vsortD(input->data, input->elementCount, order);
    }
    else {
        
        _Matrix_checkDimensionEquality(input, output);
        cblas_dcopy((UInt32)input->elementCount, input->data, 1, input->tempData, 1);
        vDSP_vsortD(input->tempData, input->elementCount, order);
        cblas_dcopy((UInt32)input->elementCount, input->tempData, 1, output->data, 1);
    }
}

OVERLOADED void Matrix_sortUsingIndices(Matrix *input, Matrix *indices)
{
    Matrix_sortUsingIndices(input, indices, input);
}

OVERLOADED void Matrix_sortUsingIndices(Matrix *input, Matrix *indices, Matrix *output)
{
    _Matrix_checkDimensionEquality(3, input, indices, output);
    
    cblas_dcopy((UInt32)input->elementCount, input->data, 1, input->tempData, 1);
    
    for (size_t i = 0; i < input->elementCount; ++i) {
        
        size_t index = indices->data[i];
        
        _Matrix_checkForValidIndex(input, index);
        
        output->data[i] = input->tempData[index];
    }
}

void detrend(Float64 *input, size_t inputLength, Float64 *output)
{
    const size_t lbp = 2;
    Float64 a1[lbp * inputLength];
    Float64 a2[lbp * inputLength];
    Float64 inputCopy1[inputLength];
    Float64 inputCopy2[inputLength];
    
    cblas_dcopy((int)inputLength, input, 1, inputCopy1, 1);
    cblas_dcopy((int)inputLength, input, 1, inputCopy2, 1);
    
    
    Float64 work[lbp * inputLength];
    
    Float64 one = 1;
    vDSP_vrampD(&one, &one, a1, 1, inputLength);
    Float64 length = inputLength;
    vDSP_vsdivD(a1, 1, &length, a1, 1, inputLength);
    vDSP_vfillD(&one, &a1[inputLength], 1, inputLength);
    cblas_dcopy((int)(inputLength), a1, 1, a2, 2);
    cblas_dcopy((int)(inputLength), &a1[inputLength], 1, &a2[1], 2);
    
    char trans = 'N';
    int n = 2;
    int m = (int)inputLength;
    int nrhs = 1;
    int lwork = (int)(lbp * inputLength);
    int info = 0;
    
    dgels_(&trans, &m, &n, &nrhs, a1, &m, inputCopy2, &m, work, &lwork, &info);
    
    vDSP_mmulD(a2, 1, inputCopy2, 1, a1, 1, m, 1, n);
    
    vDSP_vsubD(a1, 1, inputCopy1, 1, output, 1, inputLength);
}


OVERLOADED void Matrix_detrend(Matrix *input)
{
    
    detrend(input->data, input->columnCount, input->data);
}

OVERLOADED void Matrix_detrend(Matrix *input, Matrix *output)
{
    _Matrix_checkForVector(input);
    _Matrix_checkDimensionEquality(2, input, output);
    
    detrend(input->data, input->columnCount, output->data);
}

OVERLOADED void Matrix_round(Matrix *input)
{
    Matrix_round(input, input);
}

OVERLOADED void Matrix_round(Matrix *input, Matrix *output)
{
    _Matrix_checkDimensionEquality(input, output);
    
    if (input->isSubmatrix
        ||
        output->isSubmatrix) {
        
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            for (size_t j = 0; j < input->columnCount; ++j) {
                
                Matrix_getRow(output, i)[j] = round(Matrix_getRow(input, i)[j]);
            }
        }
    }
    else {
        
        for (size_t i = 0; i < input->elementCount; ++i) {
            
            output->data[i] = round(input->data[i]);
        }
    }
}

Float64 Matrix_sumIndexedElements(Matrix *input, Matrix *indices)
{
    _Matrix_checkForVector(input);
    _Matrix_checkForVector(indices);
    
    Float64 sum = 0;
    
    for (size_t i = 0; i < indices->columnCount; ++i) {
        
        size_t currentIndex = indices->data[i];
        
        _Matrix_checkForValidIndex(input, currentIndex);
        
        sum += input->data[currentIndex];
        
    }
    
    return sum;
}
OVERLOADED void Matrix_absolute(Matrix *input)
{
    Matrix_absolute(input, input);
}

OVERLOADED void Matrix_absolute(Matrix *input, Matrix *output)
{
    _Matrix_checkDimensionEquality(2, input, output);
    
    if (input->isSubmatrix
        ||
        output->isSubmatrix) {
        
        for (size_t i = 0; i < input->rowCount; ++i) {
            
            vDSP_vabsD(Matrix_getRow(input, i), 1, Matrix_getRow(output, i), 1, input->columnCount);
        }
    }
    else {
        
        vDSP_vabsD(input->data, 1, output->data, 1, input->elementCount);
    }
}

OVERLOADED void Matrix_reverse(Matrix *input)
{
    _Matrix_checkForVector(input);
    vDSP_vrvrsD(input->data, 1, input->elementCount);
}

OVERLOADED void Matrix_reverse(Matrix *input, Matrix *output)
{
    _Matrix_checkForVector(input);
    _Matrix_checkDimensionEquality(2, input, output);
    
    Matrix_copy(input, output);
    vDSP_vrvrsD(output->data, 1, output->elementCount);
}
