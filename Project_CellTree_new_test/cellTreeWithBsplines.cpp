#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <utility>
#include <algorithm>
#include <Windows.h>
#include <ctime>

using namespace std;

typedef double value_type;

#define m_degree 1;
const value_type L_limit = 0;
const value_type R_limit = 1;

template <typename T>
class BsplineNode;

template<typename T>
class CellNode;


#ifndef _NLISTS_H_
#define _NLISTS_H_

template <typename T>
class Nlists
{
public:
    vector<T> m_list;

public:
    Nlists() :m_list() {}
    Nlists(T val) :m_list(1, val) {}

    void push(T val) {
        m_list.push_back(val);
    }

    T& operator[](int i) {
        return m_list[i];
    }

    // �ҵ���ɾ��һ��Ԫ�أ������Ƿ�ɾ���ɹ�
    bool ndelete(T val) {
        bool flag = false;
        int i;
        for (i = 0; i < m_list.size();i++) {
            if (m_list[i] == val) {
                flag = true;
                break;
            }
        }
        for (; i < m_list.size() - 1; i++) {
            m_list[i] = m_list[i + 1];
        }
        if (flag) {
            m_list.pop_back();
        }
        return flag;
    }

    // ɾ����Ӧ��ŵ�Ԫ��
    bool indexDelete(int k) {
        if (k >= m_list.size()) {
            return false;
        }
        for (int i=k; i < m_list.size() - 1; i++) {
            m_list[i] = m_list[i + 1];
        }
        m_list.pop_back();
        return true;
    }

    int size() {
        return m_list.size();
    }
};

#endif



#ifndef _CELLTREE_H_
#define _CELLTREE_H_

template <typename T>
class Point
{
private:
    T _u;
    T _v;

public:
    Point(T u, T y) : _u(u), _v(y) {};
    T getU() { return _u; };
    T getV() { return _v; };
};

template <typename T>
class Partition
{
public:
    char _splitdir;         // �仯����'h'��'v'
    T _splitpoint;          // �̶�ֵ
    pair<T, T> _splitrange; // �仯��Χ
public:
    Partition(char splitdir, T splitpoint, pair<T, T> splitrange) : _splitdir(splitdir), _splitpoint(splitpoint), _splitrange(splitrange) {}
    Partition() :_splitdir(NULL), _splitpoint(0), _splitrange(0, 0) {}
    void operator =(Partition& temp) {
        _splitdir = temp._splitdir;
        _splitpoint = temp._splitpoint;
        _splitrange = temp._splitrange;
    }
    bool operator == (Partition& temp)
    {
        if (_splitdir == temp._splitdir && _splitpoint == temp._splitpoint && _splitrange == temp._splitrange)
            return true;
        else
            return false;
    }
};

#ifndef _BSPLINETREE_H_
#define _BSPLINETREE_H_

template <typename T>
class BsplineNode
{
public:
    vector<T> m_uNode;
    vector<T> m_vNode;
    Nlists<Partition<T>> m_Vsplits;
    Nlists<Partition<T>> m_Hsplits;
    Nlists<CellNode<T>*> m_cells;

public:
    BsplineNode() {}
    BsplineNode(vector<T> unode, vector<T> vnode) :m_uNode(unode), m_vNode(vnode) {}
    bool operator==(BsplineNode& temp) {
        if (m_uNode == temp.m_uNode && m_vNode == temp.m_vNode) {
            return true;
        }
        else {
            return false;
        }
    }

    void splitBF(Partition<T> split, Nlists<BsplineNode<T>*>& BSplines) {
        // ����split
        if (split._splitdir == 'h') {
            // ��¼�������ӵ�split��������а����ɾ����ע��Ӻ���ǰɾ
            vector<int> stack;
            for (int i = 0; i < m_Vsplits.size(); i++) {
                if (split._splitpoint == m_Vsplits[i]._splitpoint
                    && max(split._splitrange.first,m_Vsplits[i]._splitrange.first) 
                    <= min(split._splitrange.second,m_Vsplits[i]._splitrange.second)) {
                    split._splitrange.first = min(split._splitrange.first, m_Vsplits[i]._splitrange.first);
                    split._splitrange.second = max(split._splitrange.second, m_Vsplits[i]._splitrange.second);
                    stack.push_back(i);
                }
            }
            while (!stack.empty()) {
                if (!m_Vsplits.indexDelete(stack.back())) {
                    cout << "index delete Vsplit error!" << endl;
                    exit(2);
                }
                stack.pop_back();
            }
            m_Vsplits.push(split);
        }
        else if (split._splitdir == 'v') {
            vector<int> stack;
            for (int i = 0; i < m_Hsplits.size(); i++) {
                if (split._splitpoint == m_Hsplits[i]._splitpoint
                    && max(split._splitrange.first, m_Hsplits[i]._splitrange.first)
                    <= min(split._splitrange.second, m_Hsplits[i]._splitrange.second)) {
                    split._splitrange.first = min(split._splitrange.first, m_Hsplits[i]._splitrange.first);
                    split._splitrange.second = max(split._splitrange.second, m_Hsplits[i]._splitrange.second);
                    stack.push_back(i);
                }
            }
            while (!stack.empty()) {
                if (!m_Hsplits.indexDelete(stack.back())) {
                    cout << "index delete Hsplit error!" << endl;
                    exit(2);
                }
                stack.pop_back();
            }
            m_Hsplits.push(split);
        }
        
        //if (split._splitdir == 'h') {
        //    split._splitrange.first = max(split._splitrange.first, m_uNode[0]);
        //    split._splitrange.second = min(split._splitrange.second, m_uNode[2]);
        //    m_Vsplits.push(split);
        //}
        //else if (split._splitdir == 'v') {
        //    split._splitrange.first = max(split._splitrange.first, m_vNode[0]);
        //    split._splitrange.second = min(split._splitrange.second, m_vNode[2]);
        //    m_Hsplits.push(split);
        //}
        reduceBF(BSplines);
    }
    bool reduceBF(Nlists<BsplineNode<T>*>& BSplines) {
        vector<T> unode = m_uNode;
        vector<T> vnode = m_vNode;

        // ����ʹ�ù���split���������ɾ����ע��Ӻ���ǰɾ
        vector<int> Vstack;
        vector<int> Hstack;
        // ��vertical split���д���
        for (int i = 0; i < m_Vsplits.size(); i++) {
            if (m_Vsplits[i]._splitrange.first == unode[0] && m_Vsplits[i]._splitrange.second == unode[2] && m_Vsplits[i]._splitpoint != vnode[1]) {
                vnode.push_back(m_Vsplits[i]._splitpoint);
                Vstack.push_back(i);
            }
        }
        // ��horizontal split���д���
        for (int i = 0; i < m_Hsplits.size(); i++) {
            if (m_Hsplits[i]._splitrange.first == vnode[0] && m_Hsplits[i]._splitrange.second == vnode[2] && m_Hsplits[i]._splitpoint != unode[1]) {
                unode.push_back(m_Hsplits[i]._splitpoint);
                Hstack.push_back(i);
            }
        }

        //ɾ��ʹ�ù���split
        while (!Vstack.empty()) {
            if (!m_Vsplits.indexDelete(Vstack.back())) {
                cout << "index delete Vsplit error!" << endl;
                exit(3);
            }
            Vstack.pop_back();
        }
        while (!Hstack.empty()) {
            if (!m_Hsplits.indexDelete(Hstack.back())) {
                cout << "index delete Hsplit error!" << endl;
                exit(3);
            }
            Hstack.pop_back();
        }
        sort(unode.begin(), unode.end());
        sort(vnode.begin(), vnode.end());

        if (unode.size() > 3 || vnode.size() > 3) {
            for (int i = 0; i < unode.size() - 2; i++) {
                for (int j = 0; j < vnode.size() - 2; j++) {
                    BsplineNode<T>* bb = new BsplineNode<T>(vector<T>(unode.begin() + i, unode.begin() + i + 3), vector<T>(vnode.begin() + j, vnode.begin() + j + 3));
                    bool duplication = false;
                    // ��֧���ڵ�cell�Ž���
                    for (int i = 0; i < m_cells.size(); i++) {
                        if (m_cells[i]->_uRange.first >= bb->m_uNode[0] && m_cells[i]->_uRange.second <= bb->m_uNode[2]
                            && m_cells[i]->_vRange.first >= bb->m_vNode[0] && m_cells[i]->_vRange.second <= bb->m_vNode[2]) {
                            bb->m_cells.push(m_cells[i]);
                        }
                    }
                    // ����
                    for (int i = 0; (!duplication) && i < bb->m_cells.size(); i++) {
                        for (int j = 0; j < bb->m_cells[i]->_bsplines.size(); j++) {
                            if (*(bb->m_cells[i]->_bsplines[j]) == *bb) {
                                duplication = true;
                                break;
                            }
                        }
                    }
                    // ����ظ��˾�ɾ��
                    if (duplication) {
                        delete bb;
                    }
                    // ������Ϊһ���µĻ����������Ҽ������ϵ�ÿ��cell��
                    else {
                        BSplines.push(bb);
                        for (int i = 0; i < bb->m_cells.size(); i++) {
                            bb->m_cells[i]->_bsplines.push(bb);
                        }

                        // ��split���ƽ���
                        for (int i = 0; i < m_Vsplits.size(); i++) {
                            Partition<T> temp = m_Vsplits[i];
                            if (temp._splitpoint > bb->m_vNode[0] && temp._splitpoint < bb->m_vNode[2]
                                && temp._splitrange.first < bb->m_uNode[2] && temp._splitrange.second > bb->m_uNode[0]) {
                                temp._splitrange.first = max(bb->m_uNode[0], temp._splitrange.first);
                                temp._splitrange.second = min(bb->m_uNode[2], temp._splitrange.second);
                                bb->m_Vsplits.push(temp);
                            }
                        }
                        for (int i = 0; i < m_Hsplits.size(); i++) {
                            Partition<T> temp = m_Hsplits[i];
                            if (temp._splitpoint > bb->m_uNode[0] && temp._splitpoint < bb->m_uNode[2]
                                && temp._splitrange.first < bb->m_vNode[2] && temp._splitrange.second > bb->m_vNode[0]) {
                                temp._splitrange.first = max(bb->m_vNode[0], temp._splitrange.first);
                                temp._splitrange.second = min(bb->m_vNode[2], temp._splitrange.second);
                                bb->m_Hsplits.push(temp);
                            }
                        }
                        // ����ɹ�ϸ���ˣ�Ҫ�ͷ���ԭ���ڵ�
                        if (bb->reduceBF(BSplines)) {
                            delete bb;
                        }
                    }
                }
            }
            // ɾ��ԭ����������cell�ϵĺۼ�
            for (int i = 0; i < m_cells.size(); i++) {
                m_cells[i]->_bsplines.ndelete(this);
            }
            // ɾ��ԭ��������
            BSplines.ndelete(this);
            return true;
        }
        else {
            return false;
        }
    }


};

#endif

template <typename T>
class CellNode
{
public:
    pair<T, T> _uRange;
    pair<T, T> _vRange;

    Partition<T> _split;
    CellNode* _left;
    CellNode* _right;

public:
    Nlists<BsplineNode<T>*> _bsplines;

public:

    static int count;

    bool isCellOverlap(Partition<T>& temp)
    {
        if (temp._splitdir == 'h' && temp._splitpoint > _vRange.first && temp._splitpoint < _vRange.second)
        {
            if (temp._splitrange.first <= _uRange.second && temp._splitrange.second >= _uRange.first)
                return true;
        }
        else if (temp._splitdir == 'v' && temp._splitpoint > _uRange.first && temp._splitpoint < _uRange.second)
        {
            if (temp._splitrange.first <= _vRange.second && temp._splitrange.second >= _vRange.first)
                return true;
        }
        return false;
    }
    bool isCellSplit(Partition<T>& temp)
    {
        if (temp._splitdir == 'h' && temp._splitpoint >= _vRange.first && temp._splitpoint <= _vRange.second)
        {
            if (temp._splitrange.first <= _uRange.first && temp._splitrange.second >= _uRange.second)
                return true;
        }
        else if (temp._splitdir == 'v' && temp._splitpoint >= _uRange.first && temp._splitpoint <= _uRange.second)
        {
            if (temp._splitrange.first <= _vRange.first && temp._splitrange.second >= _vRange.second)
                return true;
        }
        return false;
    }
    bool isPointInCell(Point<T>& p)
    {
        T u = p.getU();
        T v = p.getV();
        if (u > _uRange.first && u <= _uRange.second && v >= _vRange.first && v <= _vRange.second)
            return true;
        else
            return false;
    }
    bool isLeaf()
    {
        if (_left)
            return false;
        else
            return true;
    }
    CellNode* getLeft() { return _left; };
    CellNode* getRight() { return _right; }

    CellNode* getCellNode() { return this; }

public:
    CellNode(pair<T, T> uRange, pair<T, T> vRange) : _uRange(uRange), _vRange(vRange), _left(NULL), _right(NULL), _bsplines() {}
    pair<T, T> getURange() { return _uRange; }
    pair<T, T> getVRange() { return _vRange; }

    bool splitCell(Partition<T> temp, Nlists<BsplineNode<T>*>& BSplines)
    {
        bool flag = false;

        // �����ж��ڻ����ڷ�Χ��
        if (isCellOverlap(temp)) 
        {
        // �����Ҷ�ӽڵ���ֱ�ӻ��֣�һ���ڵ�Ҫôû���ӽڵ㣬Ҫôͬʱ�������ӽڵ�
        // int a = isCellSplit(temp);
        if (isLeaf() && isCellSplit(temp))
        {
            if (temp._splitdir == 'h') {
                temp._splitrange.first = _uRange.first;
                temp._splitrange.second = _uRange.second;
            }
            else if (temp._splitdir == 'v') {
                temp._splitrange.first = _vRange.first;
                temp._splitrange.second = _vRange.second;
            }

            _split = temp;

            pair<T, T> uRange1;
            pair<T, T> uRange2;
            pair<T, T> vRange1;
            pair<T, T> vRange2;

            if (temp._splitdir == 'h')
            {
                uRange1 = _uRange;
                uRange2 = _uRange;
                vRange1.first = _vRange.first;
                vRange1.second = temp._splitpoint;
                vRange2.first = temp._splitpoint;
                vRange2.second = _vRange.second;
            }
            else if (temp._splitdir == 'v')
            {
                vRange1 = _vRange;
                vRange2 = _vRange;
                uRange1.first = _uRange.first;
                uRange1.second = temp._splitpoint;
                uRange2.first = temp._splitpoint;
                uRange2.second = _uRange.second;
            }
            _left = new CellNode<T>(uRange1, vRange1);
            _right = new CellNode<T>(uRange2, vRange2);

            // ���������ŵ������ӽڵ�Ļ������б��У��ڻ������ϳ�ȥԭ�ڵ㣬�����µ������ڵ�
            for (int i = 0; i < _bsplines.size(); i++) {
                _left->_bsplines.push(_bsplines[i]);
                _right->_bsplines.push(_bsplines[i]);

                _bsplines[i]->m_cells.ndelete(this);
                _bsplines[i]->m_cells.push(_left);
                _bsplines[i]->m_cells.push(_right);
            }

            // ��ÿ�����������л��ִ���
            for (int i = 0; i < _bsplines.size(); i++) {
                count++;
                _bsplines[i]->splitBF(temp, BSplines);
            }

            flag = true;
        }
        // �������Ҷ�ӽڵ㣬�����ĸ��ķ�Χ�ڣ�ѡ���Խ��������ӽڵ�
        else if (!isLeaf())
        {
            //if (_left->splitCell(temp, BSplines))
            //    flag = true;
            //if (_right->splitCell(temp, BSplines))
            //    flag = true;
            if (temp._splitdir != _split._splitdir) {
                if (_left->splitCell(temp, BSplines))
                    flag = true;
                if (_right->splitCell(temp, BSplines))
                    flag = true;
            }
            else if (temp._splitdir == 'v') {
                if (temp._splitpoint < _split._splitpoint) {
                    if (_left->splitCell(temp, BSplines)) {
                        flag = true;
                    }
                }
                else if (temp._splitpoint > _split._splitpoint) {
                    if (_right->splitCell(temp, BSplines)) {
                        flag = true;
                    }
                }
            }
            else {
                if (temp._splitpoint < _split._splitpoint) {
                    if (_left->splitCell(temp, BSplines)) {
                        flag = true;
                    }
                }
                else if (temp._splitpoint > _split._splitpoint) {
                    if (_right->splitCell(temp, BSplines)) {
                        flag = true;
                    }
                }
            }
        }
        
        }
        return flag;
    }

    // ͨ�������ҵ���Ӧcell
    CellNode* searchNode(Point<T>& p)
    {
        if (_left == NULL)
            return this;
        else if (_left->isPointInCell(p))
            return _left->searchNode(p);
        else
            return _right->searchNode(p);
    }

    // ��ȡҶ�ӽڵ�
    void getLeaf(vector<CellNode*>& temp)
    {
        if (isLeaf())
            temp.push_back(this);
        else
        {
            _left->getLeaf(temp);
            _right->getLeaf(temp);
        }
    }

    // ��ʾcell
    void cellShow()
    {
        cout << "cell : ( " << _uRange.first << ' ' << _uRange.second << " ) ";
        cout << "* ( " << _vRange.first << ' ' << _vRange.second << " )" << endl;
    }

    // ������cell���浽�ļ�
    void nodeSave(fstream& file)
    {
        if (isLeaf())
        {
            file << _uRange.first << ' ' << _uRange.second << ' ' << _vRange.first << ' ' << _vRange.second << endl;
        }
        else
        {
            _left->nodeSave(file);
            _right->nodeSave(file);
        }
    }

    ~CellNode()
    {
        delete _left;
        delete _right;
    }
};

template <typename T>
class CellTree : public CellNode<T>
{
public:
    CellTree(pair<T, T> uRange, pair<T, T> vRange) : CellNode<T>(uRange, vRange) {}

    // ��ȡ���ڵ�
    CellTree* getCellRoot() { return this; }

    // ���ص����ڵ�cell
    CellNode<T>* searchCell(Point<T>& p)
    {
        if (this->isPointInCell(p))
        {
            return this->searchNode(p);
        }
        else
            return NULL;
    }

    // ��ȡ����Ҷ�ӽڵ�
    vector<CellNode<T>*> getAllLeave()
    {
        vector<CellNode<T>*> result;
        this->getLeaf(result);
        return result;
    }

    void treeSave(string filename)
    {
        fstream file(filename, ios::out);
        //char projectPath[1000];
        //GetCurrentDirectory(1000, projectPath);
        if (!file)
        {
            cout << "File open error!";
            exit(1);
        }
        this->nodeSave(file);
        //cout << "Cell tree saved in " + string(projectPath) + "\\" + filename << endl;
        file.close();
    }
};

#endif


template<typename T>
class BsplineWithCell :public CellTree<T>
{
public:
    Nlists<BsplineNode<T>*> m_Bsplines;
public:
    BsplineWithCell(pair<T, T> uRange, pair<T, T> vRange, vector<T> unode, vector<T> vnode) :CellTree<T>(uRange, vRange), m_Bsplines()
    {
        BsplineNode<T>* b;
        for (int i = 0; i < unode.size() - 2; i++) {
            for (int j = 0; j < vnode.size() - 2; j++) {
                b = new BsplineNode<T>(vector<T>(unode.begin() + i, unode.begin() + i + 3), vector<T>(vnode.begin() + j, vnode.begin() + j + 3));
                b->m_cells.push(this);
                this->_bsplines.push(b);
                m_Bsplines.push(b);
            }
        }
    }

public:

    // ��Ҫ�����ӿ�
    //

    // ��ȡ����Cell
    vector<CellNode<T>*> getAllCellLeave()
    {
        return this->getAllLeave();
    }

    // ��ȡ����Bspline
    vector<BsplineNode<T>*> getAllBsplines()
    {
        vector<BsplineNode<T>*> result;
        for (int i = 0; i < m_Bsplines.size(); i++) {
            result.push_back(m_Bsplines[i]);
        }
        return result;
    }

    // �������л���������
    int countBsplinesNum()
    {
        return m_Bsplines.size();
    }

    // ��ȡһ��cell�����еĻ�����
    vector<BsplineNode<T>*> getAllBsplinesContainCell(CellNode<T>* cell)
    {
        vector<BsplineNode<T>*> result;
        for (int i = 0; i < cell->_bsplines.size(); i++) {
            result.push_back(cell->_bsplines[i]);
        }
        return result;
    }

    // ��cell���л���
    void splitCellsWithPartition(Partition<T>& partition, vector<Partition<T>>& partitions)
    {
        bool flag = this->splitCell(partition, m_Bsplines);
        while (flag) {
            flag = false;
            for (Partition<T>& p : partitions) {
                if (this->splitCell(p, m_Bsplines)) {
                    flag = true;
                }
                if (p == partition) {
                    break;
                }
            }
        }
    }



    //
};


typedef Point<value_type> PointD;

typedef Partition<value_type> PartitionD;
typedef CellTree<value_type> CellTreeD;
typedef CellNode<value_type> CellNodeD;
typedef BsplineNode<value_type> BsplineNodeD;
typedef BsplineWithCell<value_type> BsplineWithCellD;

int CellNodeD::count = 0;
//int BsplineNodeD::count = 0;


int main()
{

    //������ʼ�ռ�
    BsplineWithCellD tree({ 0,1 }, { 0,1 }, vector<value_type>{0, 0, 1, 1}, vector<value_type>{0, 0, 1, 1});

    // ������Ҫִ�еķָ�
    vector<PartitionD> divide;

    // ���뺬�зָ���Ϣ���ļ�
    fstream file("splits1.txt", ios::in);
    if (!file)
    {
        cout << "file open error!" << endl;
        exit(3);
    }
    char dir;
    value_type fixed;
    pair<value_type, value_type> range;
    int j = 0;
    while (!file.eof())
    {
        file >> dir >> fixed >> range.first >> range.second;
        file.get();
        divide.push_back(PartitionD(dir, fixed, range));
        //tree.splitBSplinesWithPartition(divide.back(), divide);
        //cout << ++j << endl;
    }
    cout << "Num of partition : " << divide.size() << endl;


    auto time_begin1 = clock();
    cout << "divide.size: " << divide.size() << endl;
    // ִ��ÿ���ָ��cell��Bspline�Ļ���
    int i = 0;
    for (PartitionD& d : divide)
    {
        //if (i == 327)
        //    break;
        /*if (i % 2 != 0)
            continue;*/
        tree.splitCellsWithPartition(d, divide);
        cout << ++i << endl;
    }
    auto time_end1 = clock();
    cout << "Split succeed!" << endl << endl;
    cout << time_end1 - time_begin1 << endl;

    cout << "cell_pass: " << CellNodeD::count << endl;
    //cout << "bspline_pass: " << BsplineNodeD::count << endl;

    auto time_begin2 = clock();
    // ��ȡ����cell
    vector<CellNodeD*> cells = tree.getAllCellLeave();
    cout << "cells.size: " << cells.size() << endl;
    // ��ÿ��cell��ȡ֧��������Ļ���������ʾ
    for (CellNodeD* cell : cells)
    {
        cout << cell->_bsplines.size();
        cout << endl;
    }

    // ��ȡ���л�����
    vector<BsplineNodeD*> bsplines = tree.getAllBsplines();
    fstream bsplinesFile("bsplinesData.txt", ios::out);
    for (BsplineNodeD* bs : bsplines) {
        bsplinesFile << bs->m_uNode[0] << ' ' << bs->m_uNode[1] << ' ' << bs->m_uNode[2] << ' ';
        bsplinesFile << bs->m_vNode[0] << ' ' << bs->m_vNode[1] << ' ' << bs->m_vNode[2] << endl;
    }
    bsplinesFile.close();
    //cout << "bsplines.size: " << bsplines.size() << endl;
    //for (BsplineNodeD* bs : bsplines)
    //{
    //    cout << bs->m_cells.size();
    //    cout << endl;
    //}
    // ��ȡ���л���������
    cout << "bsplines.size: " << tree.countBsplinesNum() << endl;

    //for (BsplineNodeD* bs : bsplines)
    //{
    //    bs->bsplineShow();
    //    cout << endl;
    //}
    auto time_end2 = clock();
    cout << time_end2 - time_begin2 << endl;



    // ��ÿ��cell��Ϣ������ļ��б���ͼ����ʾ
    tree.treeSave("saveData.txt");

    return 0;
}
