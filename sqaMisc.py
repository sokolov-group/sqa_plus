# Copyright 2009-2022 SecondQuantizationAlgebra Developers. All Rights Reserved.
#
# Licensed under the GNU General Public License v3.0;
# you may not use this file except in compliance with the License.
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# In addition, any modification or use of this software should
# cite the following paper:
#
#   E. Neuscamman, T. Yanai, and G. K.-L. Chan.
#   J. Chem. Phys. 130, 124102 (2009)
#
# Author: Eric Neuscamman <eric.neuscamman@gmail.com>
#

from sqaOptions import options
import time


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def makeTuples(n,inList):
    "Returns a list of all n-tuples that can be formed from the elements of inList.\n" + \
    "Warning: behavior may be crazy if inList consists of mutable objects."
    outList = []
    if n == 0:
        return [[]]
    if n == len(inList):
        return [inList]
    if n == 1:
        for i in range(len(inList)):
            outList.append([inList[i]])
        return outList
    tempList = []
    for i in range(len(inList)-n+1):
        tempList.append(inList[i])
    index = 0
    for i in range(len(tempList)):
        subList = makeTuples(n-1,inList[i+1:])
        for j in range(len(subList)):
            outList.append([tempList[i]])
            for k in range(len(subList[j])):
                outList[-1].append(subList[j][k])
    return outList


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def allDifferent(x):
    "Returns True if all the elements of x are different and False otherwise."
    for i in range(len(x)):
        if x[i] in x[i+1:]:
            return False
    return True


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def makePermutations(n):
    "Returns a list of all permutations of the order of the integers 0, 1, ..., n-1"
    if n == 1:
        return [[0]]
    outlist = []
    for i in range(n):
        temp = makePermutations(n-1)
        for j in range(len(temp)):
            for k in range(len(temp[j])):
                if temp[j][k] >= i:
                    temp[j][k] += 1
            outlist.append(temp[j])
            outlist[-1].insert(0,i)
    return outlist


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


def get_num_perms(ti,bi):
    """
    Returns the number of permutations between two lists of the integers 0 to N.
    """

    x = []
    for i in range(len(ti)):
        x.append([ti[i],bi[i]])
    x.sort(lambda a,b: cmp(a[1],b[1]))
    for i in range(len(x)):
        x[i] = x[i][0]
    i = 0
    n_perms = 0
    while i < len(x)-1:
        if x[i] != i:
            t = x[i]
            x[i] = x[t]
            x[t] = t
            n_perms += 1
        else:
            i += 1
    return n_perms


#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


