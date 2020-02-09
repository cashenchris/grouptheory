import group
import unittest

class WordTest(unittest.TestCase):
    def setUp(self):
        self.F=group.FGGroup(gens=['x','y'], inverses=['X','Y'])
        self.w1=group.Word([1,2,-1,-2],self.F)
        self.w2=group.Word('xyXY',self.F)
        self.w3=self.F.word('xyXY')
        self.w4=self.F.word([1,2,-1,-2])
        self.w5=self.F.word('xyyxYX')

        
    def test_equality_of_word_letters_for_equivalent_words_initialized_with_letters_and_numbers(self):
        # F=group.FGGroup(gens=['x','y'], inverses=['X','Y'])
        # w1=group.Word([1,2,-1,-2],F)
        # self.w2=group.Word('xyXY',F)
        # w3=F.word('xyXY')
        # w4=F.word([1,2,-1,-2])
        self.assertEqual(self.w1.letters,self.w2.letters)
        self.assertEqual(self.w1.letters,self.w3.letters)
        self.assertEqual(self.w1.letters,self.w4.letters)

    def test_equality_of_words_should_fail_in_base_group_class(self):
        # F=group.FGGroup(gens=['x','y'], inverses=['X','Y'])
        # w1=group.Word([1,2,-1,-2],F)
        # self.w2=group.Word('xyXY',F)
        with self.assertRaises(AssertionError):
            self.w1==self.w1
        with self.assertRaises(AssertionError):
            self.w1!=self.w1

    def test_cyclic_reduce(self):
        c = self.F.cyclic_reduce(self.w5)
        self.assertEqual(c.letters,[2,1])

    def test_cyclic_reducer(self):
        w0,w1 = self.F.cyclic_reducer(self.w5)
        self.assertEqual([w0.letters,w1.letters],[[-2,-1],[2,1]])

    # TODO test random_methods

#class SubgroupTest(unittest.TestCase):
#TODO test subgroup, FGSubgroup,

#class FPGroupTest(unittest.TestCase):
#TODO write some tests

class FreeReduceTest(unittest.TestCase):
    examples = ( ([1,2,-3],[1,2,-3]),
                 ([1,-1],[]),
                 ([4,4,2,-1,1,-2,3,-3,-4,3],[4,3]),
                 ([],[])
    )
    def test_freereduce(self):
        for unreduced,reduced in self.examples:
            self.assertEqual(group.freereduce(unreduced),reduced)


    
