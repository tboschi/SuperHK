#include <iostream>
#include <vector>

#include "TTree.h"
#include "TFile.h"

int main()
{
	TTree *t = new TTree("test", "test");
	t->SetDirectory(0);
	int a;
	t->Branch("a", &a);

	for (a = 0; a < 10; ++a)
		t->Fill();

	std::cout << t << " Loaded " << t->GetEntries() << std::endl;

	TTree *c = static_cast<TTree*>(t->Clone());
	c->SetDirectory(0);

	c->Fill();

	std::cout << c << " Clone has " << c->GetEntries() << std::endl;

	TTree *p = t->CloneTree();
	p->SetDirectory(0);

	p->Fill();
	p->Fill();

	std::cout << p << " Copy has " << p->GetEntries() << std::endl;

	delete t;

	return 0;
}
