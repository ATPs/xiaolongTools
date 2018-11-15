#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/home2/wli/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import ete3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        #fn, f_table = sys.argv[1:3]
        fn = sys.argv[1]
    except:
        print("Usage: *.py RAxML_bestTree.noGap", file=sys.stderr)
        sys.exit()

    nodes = []
    for each in ['6675', 6244, 6243, 6242]:
        nodes.append('%s_cp1' % each)
        nodes.append('%s_cp2' % each)


    new = []
    count = 0
    for line in cmn.file2lines(fn):
        count += 1
        t = ete3.Tree(line.replace('[&U]', ''))

        outgroup = t.get_common_ancestor(nodes)
        t.set_outgroup(outgroup)

        #R = t.get_midpoint_outgroup()
        #t.set_outgroup(R)
        #dist_dict = {each.name: t.get_distance(each)
        #        for each in t}
        leaves = [each for each in t]
        leaves = sorted(leaves, key = lambda x: t.get_distance(x))
        all_dist = {}
        for node in t.traverse("postorder"):
            all_dist[node] = node.dist

        for leaf in leaves:
            allnodes = []
            node = leaf
            while node:
                allnodes.append(node)
                node = node.up

            #to make the node from top to leaf
            allnodes.reverse()

            dist_left = 1
            for i, thisNode in enumerate(allnodes):
                dist = all_dist[thisNode]
                try:
                    isProcessed = thisNode.isProcessed
                    dist_left -= thisNode.dist
                except:
                    isProcessed = False

                if not isProcessed:
                    todoNodes = allnodes[i:]
                    #get each distance
                    dist_to_upper = [all_dist[node] for node in todoNodes]
                    print(thisNode.name, dist_to_upper)
                    total_dist = sum(dist_to_upper)
                    factor = dist_left / total_dist
                    print('cc:', dist_left, total_dist, factor)
                    thisNode.dist = dist * factor
                    dist_left -= thisNode.dist
                    thisNode.isProcessed = True
                    #print 'label dist node.dist, factor, dist_left'
                    print(thisNode.name, dist, thisNode.dist, factor, dist_left)


        info = t.write(format=1)
        new.append(info)

    new.append('')
    cmn.write_lines(new, fn + '.rescale2')



