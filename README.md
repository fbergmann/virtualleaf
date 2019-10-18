## Virtual Leaf
This is a fork from <https://github.com/litsol/virtualleaf> to integrate the root example code from supplement s1 of De Vos et al. (reference 2 below). 



## References
1. [VirtualLeaf: An Open-Source Framework for Cell-Based Modeling of Plant Tissue Growth and Development](http://www.plantphysiol.org/content/155/2/656)
2. [Putting Theory to the Test: Which Regulatory Mechanisms Can Drive Realistic Growth of a Root?](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003910)

## License
Since the original code of virtual leaf, and the supplement are licensed as GPL, so are the modifications made here. 

## Third Party code added
* win_getopt.h to compile with msvc with license: 

		
		 /*
		 * Copyright (c) 2002 Todd C. Miller <Todd.Miller@courtesan.com>
		 *
		 * Permission to use, copy, modify, and distribute this software for any
		 * purpose with or without fee is hereby granted, provided that the above
		 * copyright notice and this permission notice appear in all copies.
		 *
		 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
		 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
		 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
		 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
		 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
		 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
		 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
		 *
		 * Sponsored in part by the Defense Advanced Research Projects
		 * Agency (DARPA) and Air Force Research Laboratory, Air Force
		 * Materiel Command, USAF, under agreement number F39502-99-1-0512.
		 */
		/*-
		 * Copyright (c) 2000 The NetBSD Foundation, Inc.
		 * All rights reserved.
		 *
		 * This code is derived from software contributed to The NetBSD Foundation
		 * by Dieter Baron and Thomas Klausner.
		 *
		 * Redistribution and use in source and binary forms, with or without
		 * modification, are permitted provided that the following conditions
		 * are met:
		 * 1. Redistributions of source code must retain the above copyright
		 *    notice, this list of conditions and the following disclaimer.
		 * 2. Redistributions in binary form must reproduce the above copyright
		 *    notice, this list of conditions and the following disclaimer in the
		 *    documentation and/or other materials provided with the distribution.
		 *
		 * THIS SOFTWARE IS PROVIDED BY THE NETBSD FOUNDATION, INC. AND CONTRIBUTORS
		 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
		 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
		 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR CONTRIBUTORS
		 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
		 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
		 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
		 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
		 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
		 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
		 * POSSIBILITY OF SUCH DAMAGE.
		 */
